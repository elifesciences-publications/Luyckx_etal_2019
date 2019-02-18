function [data, stim, time, triggers, aborted] = Exp_trial(inEEG,randomise, practice, varargin)
% function [data, stim, time, triggers, aborted] = Exp_trial(inEEG,randomise, practice, [ppnr],[ntrials],[nblocks],[nsamp])

try
    %% DEFAULT VALUES
    
    optargs = {99 24 4 10 6};
    
    % Now put these defaults into the valuesToUse cell array,
    % and overwrite the ones specified in varargin.
    specif = find(~cellfun(@isempty,varargin)); % find position of specified arguments
    [optargs{specif}] = varargin{specif};
    
    % Place optional args in memorable variable names
    [ppnr, ntrials, nblocks, nsamp, numrange] = optargs{:};
    
    %% Initialise EEG system
    
    [portobject, portaddress,triggerlength,holdvalue,triggers] = connectEEG(inEEG,practice);
    
    %% Settings
    
    % Psychtoolbox defaults
    PsychDefaultSetup(2);
    
    % Initialize for PsychPortAudio
    InitializePsychSound(1);
    
    % Skip sync tests
    Screen('Preference', 'SkipSyncTests', 1);
    Screen('Preference', 'VisualDebugLevel', 1);
    %HideCursor;	% Hide the mouse cursor
    %ListenChar(2); % makes it so characters typed don't show up in the command window
    commandwindow;
    
    % Keys
    KbName('UnifyKeyNames');
    esc  	= KbName('Escape');
    space   = KbName('Space');
    left 	= KbName('f');
    right 	= KbName('j');
    
    alwaysKeys = ([esc, space, left, right]);
    
    RestrictKeysForKbCheck([alwaysKeys]); % restrict ppt to only use these keys
    
    %% Data path to save temporary variables
    
    datapath = cd;
    datafolder = fullfile(datapath,'/Numb_data/');
    
    %% Initialise data and setup
    
    [data, stim, time] = Exp_setup(randomise,practice,ppnr,ntrials,nblocks,nsamp,numrange);
    
    btrials = ntrials/nblocks; % number of trials per block
    
    %% Functions
    
    % Accuracy to money
    moneyBonus = @(maxBonus,cor,ntrials) length(find(cor == 1))/ntrials * maxBonus;
    
    % Transform rgb
    rgb = @(x) x/255;     % Transform RGB to 0 - 1 values
    
     %% Trial variables and open window
    
    screens         = Screen('Screens');
    screenNumber    = max(screens);
    
    % COLOUR variables
    
    % Basic colours
    col.white           = WhiteIndex(screenNumber); % Define black and white (white will be 1 and black 0). This is because
    col.black           = BlackIndex(screenNumber); % luminace values are (in general) defined between 0 and 1.
    
    % Background colours
    col.background      = col.white / 2;
    col.fix             = [.4 .4 .4].*col.white;
    
    % Category colours
    col.category(1,:)  	= rgb([230 159 0]).*col.white; % orange
    col.category(2,:)  	= rgb([86 180 233]).*col.white; % blue
    
    % //SCREEN variables
    
    % Open window and window data
    [w, windowRect]                 = PsychImaging('OpenWindow', screenNumber, col.background); % Open an on screen window and color it grey
    [screenXpixels, screenYpixels]  = Screen('WindowSize', w);     % Get the size of the on screen window in pixels
    [scr.xCenter, scr.yCenter]      = RectCenter(windowRect); % Get the centre coordinate of the window in pixels
    
    % Retrieve the maximum priority number
    topPriorityLevel = MaxPriority(w);
    Priority(topPriorityLevel);
    Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    % //TEXT VARIABLES (for instructions or other on-screen text)
    
    leftMargin      = 50;
    rightMargin     = leftMargin;
    topMargin       = 50;
    
    titleSize       = 60;
    textSize        = 26;
    stimSize        = 60;
    respSize        = 40;
    standardFont    = 'Calibri';
    
    Screen('TextFont', w, standardFont); % Font
    
    % //STIMULUS variables
    
    % Fixation circle
    scr.fixRad          = 10; % radius of the fixation circle
    
    % Frame for stimulus presentation
    stim.frameSide      = 300;
    stimRect            = [0 0 stim.frameSide stim.frameSide]; % Frame for stimulus presentation
    rectXpos            = screenXpixels * .5;
    rectYpos            = screenYpixels * .5;
    scr.rectCoord       = CenterRectOnPointd(stimRect, rectXpos, rectYpos);
    
    % Position of response cues
    dev                 = stim.frameSide/4; % deviation from center
    scr.penForcedResp   = 5;
    boxMargin           = 30;
    boxHeight           = respSize + boxMargin;
    boxWidth            = boxHeight;
    respFrame           = [0 0 boxWidth boxHeight]; % Frame for stimulus presentation
    scr.opt1Xpos        = scr.xCenter - boxWidth/2 - dev;
    scr.opt2Xpos        = scr.xCenter + boxWidth/2 + dev;
    scr.optYpos         = scr.yCenter;
    scr.lRespFrame  	= CenterRectOnPointd(respFrame, scr.opt1Xpos, scr.optYpos); % frame box
    scr.rRespFrame  	= CenterRectOnPointd(respFrame, scr.opt2Xpos, scr.optYpos); % frame box
    
    % Load feedback sounds
    positive            = psychwavread('Positive.wav');
    negative            = psychwavread('Negative.wav');
    goodTone            = [positive positive]';
    badTone             = [negative negative]';
    nrchannels          = size(goodTone,1);    
    pahandle            = PsychPortAudio('Open', [], [], 0, [], nrchannels);
    
    %% Instructions before experiment
    
    HideCursor;	% Hide the mouse cursor
    
    time.expStart = GetSecs; % start of experiment
    
    if practice == 1
        % Instructions before practice
        instructions    = 'beforepractice'; Exp_instructions;
    elseif practice == 0
        % Instructions after practice
        instructions    = 'afterpractice'; Exp_instructions;
    end
    
    %% Actual trial
    
    % Display blank (grey) screen
    Screen('TextFont', w, 'Calibri');
    Screen(w,'FillRect',col.background);    % blank screen
    Screen(w,'Flip');                   % write to screen
    WaitSecs(1);
    
    for t = 1:ntrials
        
        %% Start of a block
        
        % Instructions before each new block
        if mod(t,btrials) == 1
            instructions = 'startblock'; Exp_instructions;
        end
        
        %% Start of a trial
        
        if practice == 0
            disp(num2str(t));   % display trial number
            [triggers.time.trialPerBlock(t)] = sendTriggers(inEEG,practice,triggers.trialPerBlock+mod(t,btrials),portobject,portaddress,triggerlength,holdvalue);
        end
        
        % Adjust text
        Screen('TextStyle',w, 0); % bold
        Screen('TextSize', w , stimSize);
        
        %% Stimulus presentation
        
        % Fixation cross
        Screen('DrawDots', w, [scr.xCenter, scr.yCenter] , scr.fixRad, col.white, [], 2);
        [time.trialStart] = Screen(w,'Flip');
        [triggers.time.fixOnset(t)] = sendTriggers(inEEG,practice,triggers.fixOnset,portobject,portaddress,triggerlength,holdvalue);
        WaitSecs(time.fixDur);
        
        Screen(w,'Flip');
        [triggers.time.fixOffset(t)] = sendTriggers(inEEG,practice,triggers.fixOffset,portobject,portaddress,triggerlength,holdvalue);
        
        % Stimulus presentation
        for i = 1:nsamp
            
            % Fixation cross + number
            Screen('DrawDots', w, [scr.xCenter, scr.yCenter] , scr.fixRad, col.fix, [], 2);
            
            line1 = sprintf('%s',num2str(stim.samples(t,i)));
            rct = CenterRectOnPoint(Screen('TextBounds',w,line1),rectXpos,rectYpos);
            Screen('DrawText',w,line1,rct(1),rct(2),col.category(stim.category(t,i),:));
            
            Screen(w,'Flip');
            [triggers.time.stimOnset(t)] = sendTriggers(inEEG,practice,triggers.stimOnset+(i-1),portobject,portaddress,triggerlength,holdvalue);
            WaitSecs(time.stimDur);
            
            %ISI
            Screen('DrawDots', w, [scr.xCenter, scr.yCenter] , scr.fixRad, col.fix, [], 2);
            Screen(w,'Flip');
            [triggers.time.stimOffset(t)] = sendTriggers(inEEG,practice,triggers.stimOffset+(i-1),portobject,portaddress,triggerlength,holdvalue);
            WaitSecs(time.ISI);
        end
        
        % Response cue ISI
        Screen('DrawDots', w, [scr.xCenter, scr.yCenter] , scr.fixRad, col.white, [], 2);
        Screen(w,'Flip');
        WaitSecs(time.rISI);
        
        % Response cue
        Screen('DrawDots', w, [scr.xCenter, scr.yCenter] , scr.fixRad, col.fix, [], 2);
        respOptions(w,t,data,col,scr,respSize);
        respFrames(w,col,scr);
        [time.presEnd] = Screen(w,'Flip');
        [triggers.time.respBoxOnset(t)] = sendTriggers(inEEG,practice,triggers.respBoxOnset,portobject,portaddress,triggerlength,holdvalue);
        
        %% Response recording
        press       = 0;
        resp        = 0;
        keycode     = 0;
        
        while resp == 0
            [kdown, ~, codes] = KbCheck;  % check for key press
            
            % response deadline
            elapsed = GetSecs - time.presEnd;
            if elapsed > time.deadline
                data.r(t,1) = 0;
                break;
            end
            
            % check escape key
            if kdown==1
                if codes(esc)
                    
                    % Save data
                    if practice == 0
                        tmpname = sprintf('datatmp_%d.mat',ppnr);
                        save(fullfile(datafolder,tmpname),'data','stim','time');
                    end
                                        
                    aborted = 1;
                    sca;
                    closeStuff();
                    
                    if inEEG
                        CloseIOPort;
                    end
                    
                    return;
                end
                
                press = press+1;
                if codes(left) == 1 || codes(right) == 1
                    data.RT(t,1)        = GetSecs - time.presEnd;   % log RT
                    keycode             = find(codes==1);           % which button
                    data.keycode(t,1)   = keycode(1);               % take only first in case of simultaneous press
                    resp                = 1;
                end
            end
        end
        
        [triggers.time.response(t)] = sendTriggers(inEEG,practice,triggers.response,portobject,portaddress,triggerlength,holdvalue);
        
        %% Feedback
        
        % Feedback ISI
        Screen('DrawDots', w, [scr.xCenter, scr.yCenter] , scr.fixRad, col.fix, [], 2);
        respFrames(w,col,scr);
        chosenFrame(w,t,data,col,scr);
        Screen(w,'Flip');
        [triggers.time.fbOnset(t)] = sendTriggers(inEEG,practice,triggers.fbOnset,portobject,portaddress,triggerlength,holdvalue);
        WaitSecs(time.fbISI);
        
        % Response
        if data.keycode(t,1) == left
            data.r(t,1) = 1;
        elseif data.keycode(t,1) == right
            data.r(t,1) = 2;
        end
        
        % Correct?
        if data.xr(t,1) == data.r(t,1)
            data.cor(t,1) = 1;
        elseif data.r(t,1) == 0 % too slow
            data.cor(t,1) = -1;
        else
            data.cor(t,1) = 0;
        end
        
        % Feedback  
        if data.r(t) == -99
            Screen('TextStyle',w,1); % bold
            rct = CenterRectOnPoint(Screen('TextBounds',w,'X'),scr.xCenter,scr.yCenter);
            Screen('DrawText',w,'X',rct(1),rct(2),col.red);
        else
            Screen('DrawDots', w, [scr.xCenter, scr.yCenter], scr.fixRad, col.fix, [], 2);
        end
                
        outcomeSigns(w,t,data,col,scr,standardFont,respSize);
        respFrames(w,col,scr);
        chosenFrame(w,t,data,col,scr);
        
        if data.cor(t) == 1  
            PsychPortAudio('FillBuffer', pahandle, goodTone); % audio        
        else
            PsychPortAudio('FillBuffer', pahandle, badTone); % audio
        end
                
        [time.trialEnd] = Screen(w,'Flip');
        PsychPortAudio('Start', pahandle);
        WaitSecs(time.fbTime);
        PsychPortAudio('Stop', pahandle); % Stop playback
        
        [triggers.time.fbOffset(t)] = sendTriggers(inEEG,practice,triggers.fbOffset,portobject,portaddress,triggerlength,holdvalue);
        
        %% End of trial
                
        %data.accReward(t,1) = sum(data.corr(1:t));
        time.trialDur(t,1)  = time.trialEnd - time.trialStart;
        
        % ITI
        Screen('DrawDots', w, [scr.xCenter, scr.yCenter] , scr.fixRad, col.fix, [], 2);
        Screen(w,'Flip');
        WaitSecs(time.ITI(t));
        
        %% End of a block
        
        if mod(t,btrials) == 0
            
            WaitSecs(time.EBI);
                        
            blockacc = sum(data.block == data.block(t) & data.cor == 1);
            
            % Text on screen
            instructions = 'endblock'; Exp_instructions;
            
            if practice == 0
                % Save data from block
                tmpname = sprintf('datatmp_%d.mat',ppnr);
                save(fullfile(datafolder,tmpname),'data','stim','time');
            end
        end
        
    end
    
    %% End experiment
    
    if practice == 0
        
        % Translate points into money
        money = roundn(moneyBonus(2.5,data.cor,ntrials),-1);
        
        % Instructions at end of experiment
        instructions = 'endexp'; Exp_instructions;
        
        % Get length of experiment
        time.expEnd     = GetSecs;
        time.expDur     = time.expEnd - time.expStart;
        
        % Close screen
        sca
    end
    
    % Close all
    aborted = 0;
    closeStuff();
    
    % Close EEG port
    if inEEG
        CloseIOPort;
    end
    
catch ME
    
    warning('Something went wrong in the experiment.');
    
    tmpname = sprintf('datatmp_%d.mat',ppnr);
    save(fullfile(datafolder,tmpname),'data','stim','time');
    
    aborted = 1;
    sca;
    closeStuff();
    
    if inEEG
        CloseIOPort;
    end
    
    % Get length of experiment.
    if time.expStart > 0
        time.expEnd     = GetSecs;
        time.expDur     = time.expEnd - time.expStart;
    end
    
    rethrow(ME)
    
end

end

%% Commands for closure

function closeStuff()

ShowCursor;
ListenChar(0);
Priority(0);
PsychPortAudio('Close'); % Close the audio device

end

%% Response options

function respOptions(w,t,data,col,scr,respSize)

Screen('TextSize', w,respSize);
Screen('TextStyle',w,1); % bold

opt1 = 'O';
opt2 = 'B';

if data.rmap(t) == 1 % left = first picture
    rct = CenterRectOnPoint(Screen('TextBounds',w,opt1),scr.opt1Xpos,scr.optYpos);
    Screen('DrawText',w,opt1,rct(1),rct(2),col.category(1,:));
    rct = CenterRectOnPoint(Screen('TextBounds',w,opt2),scr.opt2Xpos,scr.optYpos);
    Screen('DrawText',w,opt2,rct(1),rct(2),col.category(2,:));
elseif data.rmap(t) == 2 % left = second picture
    rct = CenterRectOnPoint(Screen('TextBounds',w,opt1),scr.opt1Xpos,scr.optYpos);
    Screen('DrawText',w,opt2,rct(1),rct(2),col.category(2,:));
    rct = CenterRectOnPoint(Screen('TextBounds',w,opt2),scr.opt2Xpos,scr.optYpos);
    Screen('DrawText',w,opt1,rct(1),rct(2),col.category(1,:));
end

end

%% Response frames

function respFrames(w,col,scr)

Screen(w,'FrameRect',col.white,scr.lRespFrame,scr.penForcedResp);
Screen(w,'FrameRect',col.white,scr.rRespFrame,scr.penForcedResp);

end

%% Indicate chosen frame

function chosenFrame(w,t,data,col,scr)

if data.r(t,1) == 1    
    if data.rmap(t) == 1 % left = first picture
        Screen(w,'FrameRect',col.category(1,:),scr.lRespFrame,scr.penForcedResp);
    elseif data.rmap(t) == 2 % left = second picture
        Screen(w,'FrameRect',col.category(2,:),scr.lRespFrame,scr.penForcedResp);
    end
elseif data.r(t,1) == 2
    if data.rmap(t) == 1
        Screen(w,'FrameRect',col.category(2,:),scr.rRespFrame,scr.penForcedResp);
    elseif data.rmap(t) == 2
        Screen(w,'FrameRect',col.category(1,:),scr.rRespFrame,scr.penForcedResp);
    end
end

end

%% Feedback

function outcomeSigns(w,t,data,col,scr,respFont,respSize)

Screen('TextFont', w, respFont);
Screen('TextSize', w,respSize);
Screen('TextStyle',w,1); % bold

out1 = 'X';
out2 = '/';

if data.cor(t) == 1
    rct = CenterRectOnPoint(Screen('TextBounds',w,out1),scr.opt1Xpos,scr.optYpos);
    Screen('DrawText',w,'$',rct(1),rct(2),col.white);
    rct = CenterRectOnPoint(Screen('TextBounds',w,out1),scr.opt2Xpos,scr.optYpos);
    Screen('DrawText',w,'$',rct(1),rct(2),col.white);
else
    rct = CenterRectOnPoint(Screen('TextBounds',w,out2),scr.opt1Xpos,scr.optYpos);
    Screen('DrawText',w,out2,rct(1),rct(2),col.white);
    rct = CenterRectOnPoint(Screen('TextBounds',w,out2),scr.opt2Xpos,scr.optYpos);
    Screen('DrawText',w,out2,rct(1),rct(2),col.white);    
end

end

%% Connect EEG

function [portobject, portaddress,triggerlength,holdvalue,triggers] = connectEEG(inEEG,practice)

if inEEG && practice == 0
    
    IOPortfolder = 'C:\Users\EEG\Documents\MATLAB\IOPort';
    addpath(IOPortfolder);
    
    [portobject, portaddress] = OpenIOPort;
    triggerlength = 0.005; %send trigger for 5
    holdvalue     = 0;

    disp('EEG system initialised.');
    
else
    portobject      = [];
    portaddress     = [];
    triggerlength   = [];
    holdvalue       = [];
end  

% Define triggers
triggers = struct();
triggers.trialPerBlock  = 100;  % trial index
triggers.fixOnset   	= 1;    % fixation onset
triggers.fixOffset      = 2;    % fixation offset
triggers.stimOnset  	= 10;   % first donkey onset
triggers.stimOffset     = 20;   % first donkey offset
triggers.respBoxOnset   = 31;   % response boxes on screen
triggers.response       = 33;   % response
triggers.fbOnset        = 41;   % feedback onset
triggers.fbOffset       = 42;   % feedback offset

triggers.time.trialPerBlock = []; % trial index
triggers.time.fixOnset  	= [];    % fixation onset
triggers.time.fixOffset    	= [];    % fixation offset
triggers.time.stimOnset  	= [];   % first donkey onset
triggers.time.stimOffset   = [];   % first donkey offset
triggers.time.respBoxOnset  = [];   % response boxes on screen
triggers.time.response      = [];   % response
triggers.time.fbOnset       = [];   % feedback onset
triggers.time.fbOffset      = [];   % feedback offset

end

%% Send triggers

function [trigpoint] = sendTriggers(inEEG,practice,trig,portobject,portaddress,triggerlength,holdvalue)
    
if inEEG && practice == 0        
     io64( portobject, portaddress, trig); %this sends the trigger
     trigpoint = GetSecs;
     WaitSecs(triggerlength);
     io64( portobject, portaddress, holdvalue ); %this sets the trigger channel back to its hold value (0)
else
    trigpoint = 0;
end

end