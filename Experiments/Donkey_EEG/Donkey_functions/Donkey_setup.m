function [data, stim, time] = Donkey_setup(randomise,whichPhase,imgMap,varargin)
% [data, stim, time] = Donkey_setup(randomise, whichPhase, imgMap, [ppnr],[ntrials],[nblocks],[nsamp])

%% DEFAULT VALUES

optargs = {99 120 1 6};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
specif = find(~cellfun(@isempty,varargin)); % find position of specified arguments
[optargs{specif}] = varargin{specif};

% Place optional args in memorable variable names
[ppnr, ntrials, nblocks, nsamp] = optargs{:};

%% Initialise variables

[data, stim] = Donkey_initialise(ppnr,ntrials,nblocks,nsamp);

btrials = ntrials/nblocks;      % number of trials per block

makeLong = @(x) x(:); % function to squelch data

%% Get screen variables

nScreens        = Screen('Screens');
screenNumber    = max(nScreens);
screen          = Screen('Resolution',screenNumber);

%% Time variables

time.framedur       = 1/screen.hz;                  % approximate frame duration
if time.framedur == Inf
    time.framedur = .016666;
    warning('No frame duration registered. Default to 16.666 ms.');
end

time.fixDur     = .5;                   % fixation duration
time.stimDur    = .5;                   % duration of each sample
time.ISI        = .25;                  % interstimulusinterval
time.deadline   = 2;                    % response deadline
time.frameLag   = .05;                  % lag for forced choice frame
time.fbISI      = .2;                   % time between response and feedback
time.fbTime     = .5;                   % feedback time
time.ITI        = .5;                   % intertrialinterval
time.EBI        = 1;                   	% end block interval

time.expStart   = 0;                    % start of experiment
time.expEnd     = 0;                    % end of experiment
time.expDur     = 0;                    % duration of experiment
time.trialStart = 0;                    % start of trial
time.trialEnd   = 0;                    % end of trial
time.trialDur   = zeros(ntrials,1);     % duration of trial
time.presEnd    = 0;                    % end of stimulus presentation

%% Data variables

% Variables for data
data.sub            = ppnr*ones(ntrials,1);                     % subject index
data.block          = makeLong(repmat(1:nblocks,btrials,1));    % block index
data.phase          = whichPhase*ones(ntrials,1);                 % learning (1) or test (0) phase?
data.rmap           = [ones(ntrials/2,1); 2*ones(ntrials/2,1)]; % first picture with left key?
data.outcome        = zeros(ntrials,1);
data.blockReward    = zeros(ntrials,1);
data.accReward      = zeros(ntrials,1);
data.xr             = zeros(ntrials,1);

%% Stimulus variables
stim.proba      = linspace(.05,.95,nsamp);
stim.imgMap     = imgMap;
stim.allCombo   = [nchoosek(sort(imgMap),2); nchoosek(sort(imgMap,'descend'),2)];

if whichPhase == 2
    tmp                 = repmat(stim.allCombo,4,1);
    ncombinations       = length(tmp);
    stim.combo          = repmat(tmp,ntrials/ncombinations,1);
    data.forcedChoice   = repmat([zeros(ncombinations/2,1); ones(ncombinations/4,1); 2*ones(ncombinations/4,1)],ntrials/ncombinations,1);
elseif whichPhase == 3
    stim.combo          = repmat(stim.allCombo,ntrials/length(stim.allCombo),1);
end

%% Randomise everything

if randomise == 1
            
    % Response mapping 
    data.rmap           = Shuffle(data.rmap);

    % Each forced choice for each combination once
    randInd             = Shuffle([1:ntrials]);
    data.forcedChoice   = data.forcedChoice(randInd);
    stim.combo          = stim.combo(randInd,:);
    
end
   
%% Expected response for forced choice

data.xr(data.forcedChoice > 0 & data.rmap == data.forcedChoice) = 1; % left frame
data.xr(data.forcedChoice > 0 & data.rmap ~= data.forcedChoice) = 2; % right frame

end

%% Randomisation index within blocks (keep things balanced)
function randInd = randBlocks(ntrials,btrials)

randInd = zeros(1,ntrials);
for t = 1:btrials:ntrials
    randInd(t:t+btrials-1) = randperm(btrials)+(t-1);
end

end