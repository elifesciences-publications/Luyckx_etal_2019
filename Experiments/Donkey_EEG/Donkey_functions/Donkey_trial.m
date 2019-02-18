function [data, stim, time, triggers, aborted] = Donkey_trial(inEEG, randomise, whichPhase, imgMap, varargin)
% function [data, stim, time, triggers, aborted] = Donkey_trial(inEEG, randomise, whichPhase, imgMap, [ppnr],[ntrials],[nblocks],[nsamp])

try
    
    %% DEFAULT VALUES
    
    optargs = {99 45 1 6};
    
    % Now put these defaults into the valuesToUse cell array,
    % and overwrite the ones specified in varargin.
    specif = find(~cellfun(@isempty,varargin)); % find position of specified arguments
    [optargs{specif}] = varargin{specif};
    
    % Place optional args in memorable variable names
    [ppnr, ntrials, nblocks, nsamp] = optargs{:};
    
    %% Initialise EEG system
    
    [portobject, portaddress,triggerlength,holdvalue,triggers] = connectEEG(inEEG);
    
    %% Settings
    
    % Psychtoolbox defaults
    PsychDefaultSetup(2);
    
    % Initialize for PsychPortAudio
    InitializePsychSound(1);
    
    % Skip sync tests
    Screen('Preference', 'SkipSyncTests', 1);
    Screen('Preference', 'VisualDebugLevel', 1);
    HideCursor;	% Hide the mouse cursor
    ListenChar(2); % makes it so characters typed don't show up in the command window
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
    datafolder = fullfile(datapath,'/Donkey_data/');
    
    %% Initialise data and setup    
    [data, stim, time] = Donkey_setup(randomise,whichPhase,imgMap,ppnr,ntrials,nblocks,nsamp);
    
    btrials = ntrials/nblocks; % number of trials per block
    
    %% Functions
    
    % Bonus money (ADJUST)
    moneyBonus = @(maxBonus,total,singlePoints,ntrials) (total/singlePoints)/ntrials*maxBonus;
    
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
    
    % Other colours
    col.red             = [1 0 0].*col.white;
    col.forced          = rgb([204, 204, 153]).*col.white;
    col.choice          = col.white;
    col.donkeys         = rgb([[53,201,93];...
                                [174,36,146];...
                                [255,146,0];...
                                [255, 255, 0];...
                                [217,14,53];...
                                [19,178,207]]).*col.white;
    col.donkeys         = col.donkeys(stim.imgMap,:); 
    
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
    respSize        = 40;
    scoreSize       = 30;
    standardFont    = 'Calibri';
    scoreFont       = 'Century Gothic';
    
    Screen('TextFont', w, standardFont); % Font
    
    % //STIMULUS variables
    
    % Fixation circle
    scr.fixRad          = 10; % radius of the fixation circle
    
    % Background for stimulus presentation
    stim.frameSide      = 300;
    stimRect            = [0 0 stim.frameSide stim.frameSide]; % Frame for stimulus presentation
    rectXpos            = screenXpixels * .5;
    rectYpos            = screenYpixels * .5;
    scr.rectCoord       = CenterRectOnPointd(stimRect, rectXpos, rectYpos);
    
    % Forced choice frame (stimulus)
    penForcedStim       = 8;
    
    % Load images (in order of imgMap) -> probabilities are mapped on
    % different images when randomized
    imagetex = zeros(1,nsamp);
    j = 0;
    for i = stim.imgMap
        j = j+1;
        imdata = imread(['Donkey_' num2str(i)],'png'); % Load image
        imagetex(j) = Screen('MakeTexture', w, imdata); % Set up texture and rects
    end
    
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
    %PsychPortAudio('Volume', pahandle, .8); % set volume
    
    % Point indicator
    Screen('TextFont', w, scoreFont);
    Screen('TextSize', w, scoreSize);
    line1               = sprintf('%08s',num2str(0));
    scoreBounds         = Screen('TextBounds', w, line1, 0, 0);
    
    scr.boxPenWidth     = 2;
    boxMargin           = 10;
    boxHeight           = scoreSize + boxMargin;
    boxWidth            = scoreBounds(3) + boxMargin*2;
    scoreRect           = [0 0 boxWidth boxHeight]; % Frame for stimulus presentation
    scr.scoreXpos      	= screenXpixels - boxWidth/2 - rightMargin;
    scr.scoreYpos      	= textSize + 40;
    scr.moneyBoxCoord  	= CenterRectOnPointd(scoreRect, scr.scoreXpos, scr.scoreYpos); % frame box
      
    % Number of points to be won on each trial
    coins               = 13;
        
    %% Instructions before experiment
      
    % Get time of start of experiment
    time.expStart   = GetSecs; % start of experiment
    
    % Initial instruction
    instructions    = 'introduction'; Donkey_instructions;
    
    %% Actual trial
    
    % Display blank (grey) screen
    Screen(w,'FillRect',col.background);    % blank screen
    Screen(w,'Flip');                   % write to screen
    WaitSecs(1);
    
    for t = 1:ntrials
        
        %% Start of a block
        
        % Instructions before each new block
        if mod(t,btrials) == 1
            instructions = 'startblock'; Donkey_instructions;
        end
        
        %% Start of a trial
        
        % Display trial number
        disp(num2str(t));                
        [triggers.time.trialPerBlock(t)] = sendTriggers(inEEG,triggers.trialPerBlock+mod(t,btrials),portobject,portaddress,triggerlength,holdvalue);        
        
        %% Stimulus presentation
        
        % Fixation cross
        %pointsIndicator(w,col,scr,scoreFont,scoreSize,data.blockReward(t),whichPhase);
        Screen('DrawDots', w, [scr.xCenter, scr.yCenter], scr.fixRad, col.white, [], 2);
        [time.trialStart] = Screen(w,'Flip');
        [triggers.time.fixOnset(t)] = sendTriggers(inEEG,triggers.fixOnset,portobject,portaddress,triggerlength,holdvalue);
        WaitSecs(time.fixDur);
        
        %pointsIndicator(w,col,scr,scoreFont,scoreSize,data.blockReward(t),whichPhase);
        Screen(w,'Flip');
        [triggers.time.fixOffset(t)] = sendTriggers(inEEG,triggers.fixOffset,portobject,portaddress,triggerlength,holdvalue);
        WaitSecs(time.ISI);
        
        % Options
        for i = 1:2
            
            % Stimuli
            %pointsIndicator(w,col,scr,scoreFont,scoreSize,data.blockReward(t),whichPhase);
            Screen('DrawTexture', w, imagetex(stim.combo(t,i)),[],scr.rectCoord); % put image on screen
            
            % Frame for forced choice (with offset)
            if data.forcedChoice(t) == i
                Screen(w,'FrameRect',col.forced,scr.rectCoord,penForcedStim);
            end
            Screen(w,'Flip');
                if i == 1
                    [triggers.time.stim1Onset(t)] = sendTriggers(inEEG,triggers.stim1Onset,portobject,portaddress,triggerlength,holdvalue);
                elseif i == 2
                    [triggers.time.stim2Onset(t)] = sendTriggers(inEEG,triggers.stim2Onset,portobject,portaddress,triggerlength,holdvalue);                    
                end
            WaitSecs(time.stimDur);
            
            %screenshot(w,['Stimulus' num2str(i) '_learn' num2str(t) '.png']);
            
            % ISI
            %pointsIndicator(w,col,scr,scoreFont,scoreSize,data.blockReward(t),whichPhase);
            Screen(w,'Flip');
            if i == 1
                [triggers.time.stim1Offset(t)] = sendTriggers(inEEG,triggers.stim1Offset,portobject,portaddress,triggerlength,holdvalue);
            elseif i == 2
                [triggers.time.stim2Offset(t)] = sendTriggers(inEEG,triggers.stim2Offset,portobject,portaddress,triggerlength,holdvalue);
            end
            WaitSecs(time.ISI);
        end
        
        % Response screen
        respOptions(w,t,data,col,scr,respSize);        
        respFrames(w,t,data,stim,col,scr);
        %pointsIndicator(w,col,scr,scoreFont,scoreSize,data.blockReward(t),whichPhase);
        
        Screen('DrawDots', w, [scr.xCenter, scr.yCenter], scr.fixRad, col.black, [], 2);
        [time.presEnd] = Screen(w,'Flip');
        [triggers.time.respBoxOnset(t)] = sendTriggers(inEEG,triggers.respBoxOnset,portobject,portaddress,triggerlength,holdvalue);
        
        %screenshot(w,['Response_learn' num2str(t) '.png']);
        
        %% Response recording
        
        resp        = 0;
        keycode     = 0;
        
        while resp == 0
            [kdown, ~, codes] = KbCheck;  % check for key press
            
            % response deadline
            elapsed = GetSecs - time.presEnd;
            if elapsed > time.deadline
                break;
            end
            
            % check escape key
            if kdown==1
                if codes(esc)
                    
                    % Save data
                    if whichPhase == 2
                        tmpname = sprintf('datatmp_%d_learn.mat',ppnr);
                    elseif whichPhase == 3
                        tmpname = sprintf('datatmp_%d_test.mat',ppnr);
                    end
                    
                    save(fullfile(datafolder,tmpname),'data','stim','time');
                    
                    aborted = 1;
                    sca;
                    closeStuff();
                    
                    if inEEG
                        CloseIOPort;
                    end
                    
                    return;
                end
                
                if codes(left) == 1 || codes(right) == 1
                    keycode             = find(codes==1); % which button
                    data.keycode(t,1)   = keycode(1); % take only first in case of simultaneous press
                    
                    if codes(left) == 1
                        data.r(t,1)     = 1;
                    elseif codes(right) == 1
                        data.r(t,1)     = 2;
                    end
                    
                    % Check if they choose the forced choice option
                    if data.xr(t) == data.r(t) || data.forcedChoice(t) == 0
                        data.RT(t,1)        = GetSecs - time.presEnd;
                        resp = 1;
                    else
                        resp                = 0;
                        keycode             = 0;
                        data.keycode(t)     = -99;
                        data.r(t)           = -99;
                        data.RT(t)          = -99;
                    end                    
                end
            end
        end
        
        [triggers.time.response(t)] = sendTriggers(inEEG,triggers.response,portobject,portaddress,triggerlength,holdvalue);
        
        %% Feedback
                             
        % Indicate choice
        respOptions(w,t,data,col,scr,respSize);
        respFrames(w,t,data,stim,col,scr);        
        chosenFrame(w,t,data,stim,col,scr);        
        %pointsIndicator(w,col,scr,scoreFont,scoreSize,data.blockReward(t),whichPhase);
        
        Screen('DrawDots', w, [scr.xCenter, scr.yCenter], scr.fixRad, col.black, [], 2);
        Screen(w,'Flip');
        [triggers.time.fbOnset(t)] = sendTriggers(inEEG,triggers.fbOnset,portobject,portaddress,triggerlength,holdvalue);
        WaitSecs(time.fbISI);
        
        % Get probability of reward (if not too slow)
        if data.r(t) ~= -99
            data.outcome(t) = bandit(data.rmap(t),data.r(t),stim.proba(stim.combo(t,:)));
        end
        data.blockReward(t,1)   = sum(data.outcome(data.block == data.block(t))*coins);
        
        % Actual feedback
        if data.outcome(t) == 1
            PsychPortAudio('FillBuffer', pahandle, goodTone);            
        else
            PsychPortAudio('FillBuffer', pahandle, badTone);
        end
                
        outcomeSigns(w,t,data,col,scr,standardFont,respSize);
        respFrames(w,t,data,stim,col,scr);
        chosenFrame(w,t,data,stim,col,scr);
        %pointsIndicator(w,col,scr,scoreFont,scoreSize,data.blockReward(t),whichPhase);
        
        if data.r(t) == -99
            %Screen('DrawDots', w, [scr.xCenter, scr.yCenter], scr.fixRad, col.red, [], 2);
            Screen('TextStyle',w,1); % bold
            rct = CenterRectOnPoint(Screen('TextBounds',w,'X'),scr.xCenter,scr.yCenter);
            Screen('DrawText',w,'X',rct(1),rct(2),col.red);
        else
            Screen('DrawDots', w, [scr.xCenter, scr.yCenter], scr.fixRad, col.black, [], 2);
        end
        
        [time.trialEnd] = Screen(w,'Flip');
        PsychPortAudio('Start', pahandle);
        WaitSecs(time.fbTime);
        PsychPortAudio('Stop', pahandle); % Stop playback
        %screenshot(w,['Feedback_learn' num2str(t) '.png']);
        
        [triggers.time.fbOffset(t)] = sendTriggers(inEEG,triggers.fbOffset,portobject,portaddress,triggerlength,holdvalue);
        
        %% End of trial
        
        data.accReward(t,1)     = sum(data.outcome(1:t)*coins);
        time.trialDur(t,1)      = time.trialEnd - time.trialStart;
        
        % Keep score constant
        if mod(t,btrials) ~= 0
            data.blockReward(t+1) = data.blockReward(t);
        end
        
        % ITI
        Screen('DrawDots', w, [scr.xCenter, scr.yCenter], scr.fixRad, col.black, [], 2);
        %pointsIndicator(w,col,scr,scoreFont,scoreSize,data.blockReward(t),whichPhase);
        Screen(w,'Flip');
        WaitSecs(time.ITI);
        
        %% End of a block
        
        if mod(t,btrials) == 0
            
            WaitSecs(time.EBI);
                        
            % Text on screen
            instructions = 'endblock'; Donkey_instructions;
            
            % Save data
            if whichPhase == 2
                tmpname = sprintf('datatmp_%d_learn.mat',ppnr);
            elseif whichPhase == 3
                tmpname = sprintf('datatmp_%d_test.mat',ppnr);
            end
            
            save(fullfile(datafolder,tmpname),'data','stim','time');
        end
        
    end
    
    %% End experiment
    
    % Calculate money
    if whichPhase == 3
        money = moneyBonus(5,data.accReward(t),coins,ntrials);
    end
    
    % Instructions at end of experiment
    instructions = 'endexp'; Donkey_instructions;
    
    % Get length of experiment
    time.expEnd     = GetSecs;
    time.expDur     = time.expEnd - time.expStart;
    
    % Close all
    aborted = 0;
    closeStuff();
    
    if inEEG
        CloseIOPort;
    end
    
    if whichPhase == 3
        sca;
    end
    
catch ME
    
    warning('Something went wrong in the experiment.');
    
    if whichPhase == 2
        tmpname = sprintf('datatmp_%d_learn.mat',ppnr);
    elseif whichPhase == 3
        tmpname = sprintf('datatmp_%d_test.mat',ppnr);
    end
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

opt1 = 'A';
opt2 = 'B';

if data.rmap(t) == 1 % left = first picture
    rct = CenterRectOnPoint(Screen('TextBounds',w,opt1),scr.opt1Xpos,scr.optYpos);
    Screen('DrawText',w,opt1,rct(1),rct(2),col.white);
    rct = CenterRectOnPoint(Screen('TextBounds',w,opt2),scr.opt2Xpos,scr.optYpos);
    Screen('DrawText',w,opt2,rct(1),rct(2),col.white);
elseif data.rmap(t) == 2 % left = second picture
    rct = CenterRectOnPoint(Screen('TextBounds',w,opt1),scr.opt1Xpos,scr.optYpos);
    Screen('DrawText',w,opt2,rct(1),rct(2),col.white);
    rct = CenterRectOnPoint(Screen('TextBounds',w,opt2),scr.opt2Xpos,scr.optYpos);
    Screen('DrawText',w,opt1,rct(1),rct(2),col.white);
end

end

%% Response frames

function respFrames(w,t,data,stim,col,scr)

Screen(w,'FrameRect',col.white,scr.lRespFrame,scr.penForcedResp);
Screen(w,'FrameRect',col.white,scr.rRespFrame,scr.penForcedResp);

if data.xr(t) == 1
    Screen(w,'FrameRect',col.forced,scr.lRespFrame,scr.penForcedResp);
elseif data.xr(t) == 2
    Screen(w,'FrameRect',col.forced,scr.rRespFrame,scr.penForcedResp);
end

end

%% Indicate chosen frame

function chosenFrame(w,t,data,stim,col,scr)

if data.r(t,1) == 1    
    if data.rmap(t) == 1 % left = first picture
        Screen(w,'FrameRect',col.donkeys(stim.combo(t,1),:),scr.lRespFrame,scr.penForcedResp);
    elseif data.rmap(t) == 2 % left = second picture
        Screen(w,'FrameRect',col.donkeys(stim.combo(t,2),:),scr.lRespFrame,scr.penForcedResp);
    end
elseif data.r(t,1) == 2
    if data.rmap(t) == 1 % left = first picture
        Screen(w,'FrameRect',col.donkeys(stim.combo(t,2),:),scr.rRespFrame,scr.penForcedResp);
    elseif data.rmap(t) == 2 % left = second picture
        Screen(w,'FrameRect',col.donkeys(stim.combo(t,1),:),scr.rRespFrame,scr.penForcedResp);
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

if data.outcome(t) == 1
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

%% Points indicator

function pointsIndicator(w,col,scr,scoreFont,scoreSize,points,whichPhase)

if whichPhase == 3
    Screen('TextSize', w,scoreSize);
    Screen('TextStyle',w,0);
    Screen('TextFont', w, scoreFont);
    
    % Draw score
    line1 = sprintf('%08s',num2str(points));
    rct = CenterRectOnPoint(Screen('TextBounds',w,line1),scr.scoreXpos,scr.scoreYpos);
    Screen('DrawText',w,line1,rct(1),rct(2),col.white);
end

end

%% Bandit

function outcome = bandit(rmap,r,probaz)

if rmap == r
    probz = probaz(1);
else
    probz = probaz(2);
end

% Draw bandit
benchmark = rand;
if  benchmark <= probz
    outcome = 1; % points scored
else
    outcome = 0;
end

end

%% Screenshot

function screenshot(w,name)

imageArray = Screen('GetImage',w);
imwrite(imageArray, name);

end

%% Connect EEG

function [portobject, portaddress,triggerlength,holdvalue,triggers] = connectEEG(inEEG)

if inEEG
    
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
triggers.stim1Onset  	= 11;   % first donkey onset
triggers.stim1Offset    = 12;   % first donkey offset
triggers.stim2Onset  	= 21;   % second donkey onset
triggers.stim2Offset    = 22;   % second donkey offset
triggers.respBoxOnset   = 31;   % response boxes on screen
triggers.response       = 33;   % response
triggers.fbOnset        = 41;   % feedback onset
triggers.fbOffset       = 42;   % feedback offset

triggers.time.trialPerBlock = []; % trial index
triggers.time.fixOnset  	= [];    % fixation onset
triggers.time.fixOffset    	= [];    % fixation offset
triggers.time.stim1Onset  	= [];   % first donkey onset
triggers.time.stim1Offset   = [];   % first donkey offset
triggers.time.stim2Onset  	= [];   % second donkey onset
triggers.time.stim2Offset   = [];   % second donkey offset
triggers.time.respBoxOnset  = [];   % response boxes on screen
triggers.time.response      = [];   % response
triggers.time.fbOnset       = [];   % feedback onset
triggers.time.fbOffset      = [];   % feedback offset

end

%% Send triggers

function [trigpoint] = sendTriggers(inEEG,trig,portobject,portaddress,triggerlength,holdvalue)
    
if inEEG         
     io64( portobject, portaddress, trig); %this sends the trigger
     trigpoint = GetSecs;
     WaitSecs(triggerlength);
     io64( portobject, portaddress, holdvalue ); %this sets the trigger channel back to its hold value (0)
else
    trigpoint = 0;
end

end
