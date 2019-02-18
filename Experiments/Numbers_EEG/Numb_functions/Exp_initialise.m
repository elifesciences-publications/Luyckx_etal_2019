function [data, stim] = Exp_initialise(varargin)
% initialise all variables
% function [data, stim] = Exp_initialise([ppnr], [ntrials], [nblocks], [nsamp])

%% DEFAULT VALUES

optargs = {99 24 4 10 6};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
specif = find(~cellfun(@isempty,varargin)); % find position of specified arguments
[optargs{specif}] = varargin{specif};

% Place optional args in memorable variable names
[ppnr, ntrials, nblocks, nsamp, numrange] = optargs{:};

btrials = ntrials/nblocks;      % number of trials per block

%% Initialise all necessary variables

% Variables for data
data.sub        = -99*ones(ntrials,1);

data.trial      = [1:ntrials]';
data.block      = zeros(ntrials,1);
data.blocktype  = zeros(ntrials,1);
data.rmap       = zeros(ntrials,1);

data.xr         = zeros(ntrials,1);  % expected response
data.r          = -99*ones(ntrials,1); % ppt response
data.keycode    = -99*ones(ntrials,1); % code of key pressed
data.cor        = -1*ones(ntrials,1); % correct
data.RT         = -99*ones(ntrials,1); % reaction time
data.accReward  = -1*ones(ntrials,1); % accumulated reward

% Stimulus variables
stim.samples    = zeros(ntrials,nsamp);
stim.category   = zeros(ntrials,nsamp);
stim.mean       = zeros(ntrials,2);

end