function CreateRDM(s,inputfile,outputfile,do,data,params,pth)
%function CreateRDM(s,inputfile,outputfile,do,data,params,pth)
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

%% Indices

idx     = data.sub == params.submat(s);
cortmp  = data.RT(idx) > 0;
bindx   = cortmp & goodtrials'; % behavioural index

eindx   = cortmp(logical(goodtrials')) == 1; % eeg index
eindx   = makeLong(repmat(eindx,1,params.nsamp)');

%% Resize data

tmpsamp     = data.conds(idx,:); % conditions for RDM
condvec     = makeLong(tmpsamp(bindx,:)');

% Reshape eeg data
eegdat      = zscore(eeg.data(:,:,eindx),[],3);

%% Run regression

fprintf('\nCalculating betas subject %d.\n',params.submat(s));

results = rdm_computeConditionCoeffs(eegdat,condvec);

%% Compute RDM

fprintf('\nCalculating RDM for subject %d.\n',params.submat(s));

rdmSet  = rdm_computeRDM(results,params.disttype);
RDM     = permute(rdmSet,[2 3 1]);

%% Save RDM

if do.saveRDM
    
    rdm.data        = RDM;
    rdm.timepoints  = eeg.timepoints;
    rdm.conds       = {'1','2','3','4','5','6'};
    
    save(fullfile(pth.RDM,outputfile),'rdm');
    fprintf('\nRDMs saved for subject %d.\n',params.submat(s));
end

end

