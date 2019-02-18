function [data, stim, time, triggers, aborted] = Donkey_click(inEEG,randomise, imgMap, varargin)
% function [data, stim, time, triggers, aborted] = Donkey_click(inEEG,randomise, imgMap, [ppnr],[ntrials],[nblocks],[nsamp], [whichPhase])

try
    
    %% DEFAULT VALUES
    
    optargs = {99 60 1 6 1};
    
    % Now put these defaults into the valuesToUse cell array,
    % and overwrite the ones specified in varargin.
    specif = find(~cellfun(@isempty,varargin)); % find position of specified arguments
    [optargs{specif}] = varargin{specif};
    
    % Place optional args in memorable variable names
    [ppnr, ntrials, nblocks, nsamp, whichPhase] = optargs{:};
    
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
    
    alwaysKeys = ([esc, space]);
    
    RestrictKeysForKbCheck([alwaysKeys]); % restrict ppt to only use these keys
    
    %% Data path to save temporary variables    
    datapath = cd;
    datafolder = fullfile(datapath,'/Donkey_data/');
    
    %% Initialise data and setup    
    [data, stim, time] = Donkey_clickSetup(randomise,whichPhase,imgMap,ppnr,ntrials,nblocks,nsamp);
    
    btrials = ntrials/nblocks; % number of trials per block
    
    %% Functions
        
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
    col.choice          = col.white;
    col.donkeys         = rgb([...
        [53,201,93];...
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
    
    % Background for stimulus presentation
    stim.frameSide      = 250;
    stimRect            = [0 0 stim.frameSide stim.frameSide]; % Frame for stimulus presentation
    rectXpos            = screenXpixels * .5;
    rectYpos            = screenYpixels * .5;
    scr.rectCoord       = CenterRectOnPointd(stimRect, rectXpos, rectYpos);
    
    % Load images (in order of imgMap) -> probabilities are mapped on
    % different images when randomized
    imagetex = zeros(1,nsamp);
    j = 0;
    for i = stim.imgMap
        j = j+1;
        imdata = imread(['Donkey_' num2str(i)],'png'); % Load image
        imagetex(j) = Screen('MakeTexture', w, imdata); % Set up texture and rects
    end
    
    % Determine positions of images
    interspace      = 50;
    Xdisposition    = repmat([-1*stim.frameSide - interspace 0 stim.frameSide + interspace],1,2);
    Ydisposition    = [(-1*stim.frameSide/2 - interspace)*ones(1,3) (stim.frameSide/2 + interspace)*ones(1,3)];
        
    for i = 1:nsamp
        rectXpos = screenXpixels * .5 + Xdisposition(i);
        rectYpos = screenYpixels * .5 + Ydisposition(i);
        scr.rectCoord(:,:,:,:,i) = CenterRectOnPointd(stimRect, rectXpos, rectYpos);
    end
    
    scr.penChosen   = 6; % frame width chosen option
    
    % Mask
    col.mask        = [.5 .5 .5 .65];
    
%     % Position of tallies
%     for i = 1:nsamp
%         coords  = scr.rectCoord(:,:,:,:,i);
%         xPos    = (coords(1)+coords(3))/2;
%         yPos    = ((coords(2)+coords(4)) + interspace + stim.frameSide)/2;
%         
%         scr.tallyXcoords(i) = xPos;
%         scr.tallyYcoords(i) = yPos;
%     end
    
    % Load feedback sounds
    positive            = psychwavread('Positive.wav');
    negative            = psychwavread('Negative.wav');
    goodTone            = [positive positive]';
    badTone             = [negative negative]';
    nrchannels          = size(goodTone,1);    
    pahandle            = PsychPortAudio('Open', [], [], 0, [], nrchannels);
    %PsychPortAudio('Volume', pahandle, .8); % set volume
      
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
    ShowCursor;
    
    for t = 1:ntrials
        
        %% Start of a block
        
        % Instructions before each new block
        if mod(t,btrials) == 1
            instructions = 'startblock'; Donkey_instructions;
        end
        
        %% Start of a trial
        
        % Display trial number
        disp(num2str(t));
        
        % Timing start trial
        time.trialStart = GetSecs;
        
        %% Stimulus presentation
        
        % All stimuli
        for i = 1:nsamp
            Screen('DrawTexture', w, imagetex(stim.locations(i)),[],scr.rectCoord(:,:,:,:,i)); % put image on screen
        end
        
        % Mask for depleted ones
        if any(data.tally(t,:) == 0)
            masked = find(ismember(stim.locations,find(data.tally(t,:) == 0)));
            for i = masked
                Screen('FillRect', w, col.mask,scr.rectCoord(:,:,:,:,i)); % put mask
            end
        end
        
        Screen(w,'Flip');  
        [triggers.time.respBoxOnset(t)] = sendTriggers(inEEG,triggers.respBoxOnset,portobject,portaddress,triggerlength,holdvalue);
        %screenshot(w,'Stimuli_click.png');
        
        %% Response recording
        
        resp        = 0;
        
        while resp == 0
            
            % Check escape key
            [kdown, ~, codes] = KbCheck;  % check for key press
            
            if kdown==1
                if codes(esc)
                    
                    % Save data
                    tmpname = sprintf('datatmp_%d_click.mat',ppnr);
                    save(fullfile(datafolder,tmpname),'data','stim','time');
                    
                    aborted = 1;
                    sca;
                    closeStuff();
                    
                    if inEEG
                        CloseIOPort;
                    end
                    
                    return;
                end
            end
            
            % Track mouse
            inside          = logical(zeros(1,nsamp));
            [x, y, buttons] = GetMouse(w);
            for i = 1:nsamp
                inside(i) = IsInRect(x, y, scr.rectCoord(:,:,:,:,i));
            end
            
            if any(inside) && buttons(1) % check if click was inside a box
                resp                = 1;
                data.chosenLocs(t)  = find(inside == true); % chosen location
                data.chosen(t)      = stim.locations(inside == true); % probability chosen
                data.RT(t)          = GetSecs - time.trialStart;
                
                if data.tally(t,data.chosen(t)) == 0
                    resp        = 0;
                    data.RT(t)  = -99;
                end               
            end 
        end
              
        [triggers.time.response(t)] = sendTriggers(inEEG,triggers.response,portobject,portaddress,triggerlength,holdvalue);
        
        %% Feedback
              
        WaitSecs(time.fbISI);
        
        % Get probability of reward (if not too slow)
        data.outcome(t) = stim.bandit(data.tally(t,data.chosen(t)),data.chosen(t));
              
        % All stimuli
        allOpts     = 1:6;
        leftOvers   = allOpts(~ismember(1:nsamp,data.chosenLocs(t)));
        for i = leftOvers
            Screen('DrawTexture', w, imagetex(stim.locations(i)),[],scr.rectCoord(:,:,:,:,i)); % put images on screen
        end
        
        % Mask for depleted ones
        if any(data.tally(t,:) == 0)
            masked = find(ismember(stim.locations,find(data.tally(t,:) == 0)));
            for i = masked
                Screen('FillRect', w, col.mask,scr.rectCoord(:,:,:,:,i)); % put mask
            end
        end
        
        Screen(w,'FrameRect',col.donkeys(data.chosen(t),:),scr.rectCoord(:,:,:,:,data.chosenLocs(t)),scr.penChosen);
        %tallyIndicator(w,t,data,stim,col,scr,scoreFont,scoreSize);
        
        Screen(w,'Flip');
        [triggers.time.fbOnset(t)] = sendTriggers(inEEG,triggers.fbOnset,portobject,portaddress,triggerlength,holdvalue);
        WaitSecs(time.fbISI);

        % Actual feedback
        if data.outcome(t) == 1
            PsychPortAudio('FillBuffer', pahandle, goodTone);            
        else
            PsychPortAudio('FillBuffer', pahandle, badTone);
        end
           
        % All stimuli
        for i = leftOvers
            Screen('DrawTexture', w, imagetex(stim.locations(i)),[],scr.rectCoord(:,:,:,:,i)); % put images on screen
        end
        
        % Mask for depleted ones
        if any(data.tally(t,:) == 0)
            masked = find(ismember(stim.locations,find(data.tally(t,:) == 0)));
            for i = masked
                Screen('FillRect', w, col.mask,scr.rectCoord(:,:,:,:,i)); % put mask
            end
        end
        
        Screen(w,'FrameRect',col.donkeys(data.chosen(t),:),scr.rectCoord(:,:,:,:,data.chosenLocs(t)),scr.penChosen);
        
        outcomeSigns(w,t,data,col,scr,standardFont,respSize);
        %tallyIndicator(w,t,data,stim,col,scr,scoreFont,scoreSize);
        
        [time.trialEnd] = Screen(w,'Flip');
        [triggers.time.fbSignOnset(t)] = sendTriggers(inEEG,triggers.fbSignOnset,portobject,portaddress,triggerlength,holdvalue);
        %screenshot(w,'Feedback_click.png');
        
        PsychPortAudio('Start', pahandle);
        WaitSecs(time.fbTime);
        PsychPortAudio('Stop', pahandle); % Stop playback
        
        [triggers.time.fbOffset(t)] = sendTriggers(inEEG,triggers.fbOffset,portobject,portaddress,triggerlength,holdvalue);
        
        %% End of trial
        
        % Adjust tally
        if t < ntrials
            data.tally(t+1,:)               = data.tally(t,:);
            data.tally(t+1,data.chosen(t))  = data.tally(t,data.chosen(t))-1;
        end
        
        % Timing trial
        time.trialDur(t,1)      = time.trialEnd - time.trialStart;
                
        % ITI
        WaitSecs(time.ITI);
        
    end
    
    %% End experiment
    
    WaitSecs(time.EBI);
    
    % Save temporary file before things crash at the end
    tmpname = sprintf('datatmp_%d_click.mat',ppnr);
    save(fullfile(datafolder,tmpname),'data','stim','time');
        
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
    
    %sca;
    
catch ME
    
    warning('Something went wrong in the experiment.');
    
    tmpname = sprintf('datatmp_%d_click.mat',ppnr);
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

%% Feedback

function outcomeSigns(w,t,data,col,scr,respFont,respSize)

Screen('TextFont', w, respFont);
Screen('TextSize', w,respSize);
Screen('TextStyle',w,1); % bold

out1 = 'X';
out2 = '/';

% Coords
coords  = scr.rectCoord(:,:,:,:,data.chosenLocs(t));
xPos    = (coords(1)+coords(3))/2;
yPos    = (coords(2)+coords(4))/2;

if data.outcome(t) == 1
    rct = CenterRectOnPoint(Screen('TextBounds',w,out1),xPos,yPos);
    Screen('DrawText',w,'$',rct(1),rct(2),col.white);
else
    rct = CenterRectOnPoint(Screen('TextBounds',w,out2),xPos,yPos);
    Screen('DrawText',w,out2,rct(1),rct(2),col.white);   
end

end

%% Indicate running tally and lock

function tallyIndicator(w,t,data,stim,col,scr,scoreFont,scoreSize)

Screen('TextFont', w,scoreFont);
Screen('TextSize', w,scoreSize);
Screen('TextStyle',w,1);

for i = 1:size(data.tally,2)   
    
    line1 = num2str(data.tally(t,i));    
    rct = CenterRectOnPoint(Screen('TextBounds',w,line1),scr.tallyXcoords(stim.locations == i),scr.tallyYcoords(stim.locations == i));
    Screen('DrawText',w,line1,rct(1),rct(2),col.white);    
    
    if (data.locked(t,i) > 0 || data.stageLock(t,i) == 1) && data.tally(t,i) ~=0
        Screen('DrawText',w,'X',rct(1),rct(2),col.white);
    end
    
end

end

%% Bandit

function outcome = bandit(chosen,probaz)

% Draw bandit
benchmark = rand;
if  benchmark <= probaz(chosen)
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
triggers.respBoxOnset   = 31;   % response boxes on screen
triggers.response       = 33;   % response
triggers.fbOnset        = 41;   % feedback onset
triggers.fbOffset       = 42;   % feedback offset
triggers.fbSignOnset    = 45;   % onset of sign in frame

triggers.time.trialPerBlock = []; % trial index
triggers.time.fixOnset  	= [];    % fixation onset
triggers.time.fixOffset    	= [];    % fixation offset
triggers.time.respBoxOnset  = [];   % response boxes on screen
triggers.time.response      = [];   % response
triggers.time.fbOnset       = [];   % feedback onset
triggers.time.fbOffset      = [];   % feedback offset
triggers.time.fbSignOnset   = [];   % onset of sign in frame


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

