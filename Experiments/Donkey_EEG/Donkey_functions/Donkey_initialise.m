function [data, stim] = Donkey_initialise(varargin)
% initialise all variables
% function [data, stim] = Donkey_initialise([ppnr], [ntrials], [nblocks], [nsamp])

%% DEFAULT VALUES
optargs = {99 60 1 6};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
specif = find(~cellfun(@isempty,varargin)); % find position of specified arguments
[optargs{specif}] = varargin{specif};

% Place optional args in memorable variable names
[ppnr, ntrials, nblocks, nsamp] = optargs{:};

btrials = ntrials/nblocks;      % number of trials per block

%% Initialise all necessary variables

% Variables for data
data.sub        = -99*ones(ntrials,1);

data.trial      = [1:ntrials]';
data.block      = zeros(ntrials,1);
data.phase      = -99*ones(ntrials,1);

data.forcedChoice = zeros(ntrials,1);

data.outcome    = -99*ones(ntrials,1);
data.rmap       = zeros(ntrials,1); % which response on left?
data.r          = -99*ones(ntrials,1); % ppt response
data.RT         = -99*ones(ntrials,1); % reaction time
data.xr         = -99*ones(ntrials,1); % expected response when forced
data.keycode    = -99*ones(ntrials,1); % code of key pressed
data.RT         = -99*ones(ntrials,1); % reaction time
data.blockReward = -99*ones(ntrials,1); % reward per block
data.accReward  = -99*ones(ntrials,1); % accumulated reward

% Stimulus variables
stim.proba      = zeros(1,nsamp);
stim.imgMap     = zeros(1,nsamp);
stim.allCombo   = zeros(nsamp^2-nsamp,2); % all possible combinations of pairs
stim.combo      = zeros(ntrials,2); % combination of bandit on each trial
stim.frameSide  = 0; % length/width of frame for stimulus presentation

end