function CreateRDM_bandits(s,inputfile,outputfile,do,data,params,pth)
%function CreateRDM_bandits(s,inputfile,outputfile,do,data,params,pth)
%   Detailed explanation goes here

%% Initialise

% Index for exluding trials
goodtrials  = [];
bindx       = []; % behavioural index
eindx       = []; % eeg index

% Functions
makeLong = @(x) x(:);

%% Load data

% Load EEG data
load(fullfile(pth.raw,sprintf('%sEEG_samples.mat',inputfile)));

% Load bad trials
load(fullfile(pth.raw,sprintf('%srejectedTrials.mat',inputfile)));

%% Good trials

goodtrials  = 1-rejectedTrialz(1:params.ntrials);

% Variables
ntimepoints = length(eeg.timepoints); % number of time points in epoch
rtrials     = length(find(goodtrials == 1)); % number of preserved trials

%% Iterate over bandits

for b = 1:2
    
    %% Indices
    
    idx     = data.sub == params.submat(s);
    cortmp  = data.RT(idx) > 0;
    bindx   = cortmp & goodtrials'; % behavioural index
    
    eindx   = cortmp(logical(goodtrials')) == 1; % eeg index
    eindx   = makeLong(repmat(eindx,1,params.nsamp)');
    
    if b == 1
        eindx(2:2:length(eindx)) = 0; % take out even epochs (second bandit)
    elseif b == 2
        eindx(1:2:length(eindx)) = 0; % take out odd epochs (first bandit)
    end
    
    %% Resize data
    
    tmpsamp     = data.conds(idx,b); % perceived ranks
    condvec     = makeLong(tmpsamp(bindx,:)');
    
    % Reshape eeg data
    eegdat      = zscore(eeg.data(:,:,eindx),[],3);
    
    %% Run regression
    
    fprintf('\nCalculating betas subject %d, bandit %d.\n',params.submat(s),b);
    
    results = rdm_computeConditionCoeffs(eegdat,condvec);
    
    %% Compute RDM
    
    fprintf('\nCalculating RDM for subject %d.\n',params.submat(s));
    
    rdmSet  = rdm_computeRDM(results,params.disttype);
    RDM     = permute(rdmSet,[2 3 1]);
    
    %% Save RDM
    
    if do.saveRDM
        
        rdm.data        = RDM;
        rdm.timepoints  = eeg.timepoints;
        rdm.conds       = {'b1','b2','b3','b4','b5','b6'};
        
        save(fullfile(pth.RDM,sprintf('%s_bandit%d',outputfile,b)),'rdm');
        fprintf('\nRDMs saved for subject %d.\n',params.submat(s));
    end
    
end

end

