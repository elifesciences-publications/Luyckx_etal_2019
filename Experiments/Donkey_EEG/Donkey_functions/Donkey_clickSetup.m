function [data, stim, time] = Donkey_clickSetup(randomise,whichPhase,imgMap,varargin)
% [data, stim, time] = Donkey_setup(randomise, whichPhase, imgMap, [ppnr],[ntrials],[nblocks],[nsamp])

%% DEFAULT VALUES

optargs = {99 30 1 6};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
specif = find(~cellfun(@isempty,varargin)); % find position of specified arguments
[optargs{specif}] = varargin{specif};

% Place optional args in memorable variable names
[ppnr, ntrials, nblocks, nsamp] = optargs{:};

%% Initialise variables

btrials = ntrials/nblocks;      % number of trials per block
nclicks = ntrials/nsamp;        % number of clicks possible for each bandit

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

time.fbISI      = .2;                  % time between response and feedback
time.fbTime     = .5;                   % feedback time
time.ITI        = .1;                   % intertrialinterval
time.EBI        = 1;                   	% end block interval

time.expStart   = 0;                    % start of experiment
time.expEnd     = 0;                    % end of experiment
time.expDur     = 0;                    % duration of experiment
time.trialStart = 0;                    % start of trial
time.trialEnd   = 0;                    % end of trial
time.trialDur   = zeros(ntrials,1);     % duration of trial

%% Data variables

% Variables for data
data.sub            = ppnr*ones(ntrials,1);                     % subject index
data.block          = makeLong(repmat(1:nblocks,btrials,1));    % block index
data.phase          = whichPhase*ones(ntrials,1);
data.RT             = -99*ones(ntrials,1);
data.chosen         = zeros(ntrials,1);
data.chosenLocs     = zeros(ntrials,1);
data.outcome        = zeros(ntrials,1);
data.tally          = ntrials/nsamp*ones(ntrials,nsamp);

%% Stimulus variables
stim.proba          = linspace(.05,.95,nsamp);
stim.imgMap         = imgMap;
stim.locations      = 1:nsamp; % location of images, numbers indicate which probability is on position 1, etc ... example: loc(1) = 4 => 4th probability is on position 1
stim.bandit         = zeros(nclicks,nsamp);

%% Randomise everything

if randomise == 1
 
    % Randomise locations of images on screen
    stim.locations = Shuffle(stim.locations);
    
    % Outcome of donkeys
    for s = 1:nsamp
        winz = round(stim.proba(s) * nclicks);
        stim.bandit(:,s) = Shuffle([ones(winz,1); zeros(nclicks-winz,1)]);
    end
    
end

end