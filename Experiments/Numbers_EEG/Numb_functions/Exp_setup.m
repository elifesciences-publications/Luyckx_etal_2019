function [data, stim, time] = Exp_setup(randomise, varargin)
% [data, stim, time] = Exp_setup(randomise, [practice], [ppnr],[ntrials],[nblocks],[nsamp])

%% DEFAULT VALUES

optargs = {0 99 24 4 10 6};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
specif = find(~cellfun(@isempty,varargin)); % find position of specified arguments
[optargs{specif}] = varargin{specif};

% Place optional args in memorable variable names
[practice, ppnr, ntrials, nblocks, nsamp, numrange] = optargs{:};

%% Initialise variables

[data, stim] = Exp_initialise(ppnr,ntrials,nblocks,nsamp,numrange);

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
    warning('No frame duration registered.');
end

time.fixDur     = .5;                   % fixation duration
time.stimDur    = 17*time.framedur; 	% duration of each sample
time.ISI        = 4*time.framedur;      % interstimulusinterval
time.rISI       = time.stimDur + time.ISI; % response cue interval
time.fbISI      = .1;                   % before feedback
time.fbTime     = time.rISI;            % feedback time
time.ITI        = .5 + rand(ntrials,1); % intertrialinterval
time.EBI        = 1;                   	% end block interval

time.expStart   = 0;                    % start of experiment
time.expEnd     = 0;                    % end of experiment
time.expDur     = 0;                    % duration of experiment
time.trialStart = 0;                    % start of trial
time.trialEnd   = 0;                    % end of trial
time.trialDur   = zeros(ntrials,1);     % duration of trial
time.presEnd    = 0;                    % end of sample presentation
time.deadline   = 2;                    % response deadline

%% Data variables

% Variables for data
data.sub    = ppnr*ones(ntrials,1);                     % subject index
data.block  = makeLong(repmat(1:nblocks,btrials,1));    % block index

if rem(ppnr,2) == 1
    data.rmap = ones(ntrials,1); % left = cat1, right = cat2
elseif rem(ppnr,2) == 0
    data.rmap = 2*ones(ntrials,1); % left = cat2, right = cat1
end

%% Stimulus variables

% Create stimulus data
stim.samples    = randi([1 numrange],ntrials,nsamp);
stim.category   = [ones(ntrials,nsamp/2) 2*ones(ntrials,nsamp/2)];

%% Randomise everything

if randomise == 1
        
    % Generate trials      
    for t = 1:ntrials
        
        meanz   = zeros(1,2);        
        i       = 1;
        
        while i < 1000 && meanz(1) == meanz(2)
            samplez     = randi([1 numrange],1,nsamp);
            catz        = Shuffle(stim.category(1,:));
            
            meanz(1)    = mean(samplez(catz == 1));
            meanz(2)    = mean(samplez(catz == 2));
            i = i+1;
        end   
        
        stim.samples(t,:)  	= samplez;
        stim.category(t,:) 	= catz;
        stim.mean(t,:)      = meanz;
        
    end
    
else
    for t = 1:ntrials
        stim.mean(t,1) = mean(stim.samples(t,stim.category(t,:) == 1));
        stim.mean(t,2) = mean(stim.samples(t,stim.category(t,:) == 2));
    end
end
 
%% Calculate rest

data.xr(stim.mean(:,1) > stim.mean(:,2) & data.rmap == 1) = 1;
data.xr(stim.mean(:,1) < stim.mean(:,2) & data.rmap == 1) = 2;
data.xr(stim.mean(:,1) > stim.mean(:,2) & data.rmap == 2) = 2;
data.xr(stim.mean(:,1) < stim.mean(:,2) & data.rmap == 2) = 1;

end