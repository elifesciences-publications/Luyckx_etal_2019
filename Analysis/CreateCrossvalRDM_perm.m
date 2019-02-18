function CreateCrossvalRDM_perm(combo,num_data,donk_data,do,params,num_paths,donk_paths)
%function CreateCrossvalRDM_perm(combo,num_data,donk_data,do,params,num_paths,donk_paths)

s = combo(1);
z = combo(2);

subvec  = setdiff(1:params.nsubj,s);
submat2 = params.submat(subvec);

%% Initialise Numbers

% Index for exluding trials
goodtrials  = [];
bindx       = []; % behavioural index
eindx       = []; % eeg index

% Functions
makeLong = @(x) x(:);

% Load data
inputfile   = sprintf('Numbers_%03d_',params.submat(s));

% Load EEG data
load(fullfile(num_paths.data.EEG.raw,sprintf('%sEEG_samples.mat',inputfile)));

% Load bad trials
load(fullfile(num_paths.data.EEG.raw,sprintf('%srejectedTrials.mat',inputfile)));

% Good trials
goodtrials  = 1-rejectedTrialz;

% Variables
timepoints1 = eeg.timepoints;
ntimepoints = length(timepoints1); % number of time points in epoch
rtrials     = length(find(goodtrials == 1)); % number of preserved trials

% Indices
idx     = num_data.sub == params.submat(s);
cortmp  = num_data.cor(idx) >= 0;
bindx   = cortmp & goodtrials'; % behavioural index

eindx   = cortmp(logical(goodtrials')) == 1; % eeg index
eindx   = makeLong(repmat(eindx,1,size(params.num_conds,2))'); % expand from n trials to n epochs

% Resize data

% Create conditions vector
tmpsamp     = params.num_conds(idx,:);
condvec1    = makeLong(tmpsamp(bindx,:)');

% Reshape eeg data
eegdat1     = eeg.data(:,:,eindx);

%% Initialise Donkey

% Index for exluding trials
goodtrials  = [];
bindx       = []; % behavioural index
eindx       = []; % eeg index

% Load data Donkey
inputfile   = sprintf('Donkey_%03d_',submat2(z));

% Load EEG data
load(fullfile(donk_paths.data.EEG.raw,sprintf('%s%s_EEG_samples.mat',inputfile,params.whichPhase)));

% Load bad trials
load(fullfile(donk_paths.data.EEG.raw,sprintf('%s%s_rejectedTrials.mat',inputfile,params.whichPhase)));

% Good trials Donkey
goodtrials  = 1-rejectedTrialz;

% Variables
timepoints2 = eeg.timepoints;
ntimepoints = length(timepoints2); % number of time points in epoch
rtrials     = length(find(goodtrials == 1)); % number of preserved trials

% Indices Donkey
idx     = donk_data.sub == submat2(z);
cortmp  = donk_data.RT(idx) > 0;
bindx   = cortmp & goodtrials'; % behavioural index

eindx   = cortmp(logical(goodtrials')) == 1; % eeg index
eindx   = makeLong(repmat(eindx,1,size(params.donk_conds,2))'); % expand from n trials to n epochs

% Resize data Donkey
tmpsamp 	= params.donk_conds(idx,:); % perceived ranks RL
condvec2  	= makeLong(tmpsamp(bindx,:)');

% Reshape eeg data
eegdat2     = eeg.data(:,:,eindx);

%% Concatenate two datasets

mutual_time = intersect(timepoints1,timepoints2);
eeg1wdw     = timepoints1 >= mutual_time(1) & timepoints1 <= mutual_time(end);
eeg2wdw     = timepoints2 >= mutual_time(1) & timepoints2 <= mutual_time(end);

numb_sz     = size(eegdat1(:,eeg1wdw,:));
donk_sz     = size(eegdat2(:,eeg2wdw,:));

numbeeg     = zscore(eegdat1(:,eeg1wdw,:),[],3); % z-score data
donkeeg     = zscore(eegdat2(:,eeg2wdw,:),[],3); % z-score data

%% Compute betas

fprintf('\nCalculating betas subject %d.\n',params.submat(s));

% Run regression for numbers
results_num  = rdm_computeConditionCoeffs(numbeeg,condvec1);

% Run regression Donkeys
results_donk = rdm_computeConditionCoeffs(donkeeg,condvec2);

%% Compute RDM

fprintf('\nCalculating cross-validation RDM for subject %d to sub %d.\n',params.submat(s),submat2(z));

results.betas   = cat(3,results_num.betas,results_donk.betas);
results.resids  = cat(3,results_num.resids,results_donk.resids);

RDM             = rdm_computeRDM_crossval(results,params.disttype);

%% Save RDM

if do.saveRDM
    
    rdm.data        = RDM;
    rdm.timepoints  = mutual_time;
    rdm.conds       = {'1','2','3','4','5','6','b1','b2','b3','b4','b5','b6'};
    
    outputfile  = sprintf('Crossval_%03d_RDM_%s_%s_perm%d',params.submat(s),params.whichPhase,submat2(z),params.disttype,submat2(z));
    
    save(fullfile(donk_paths.data.save,outputfile),'rdm');
    fprintf('\nRDMs saved for subject %d to sub %d.\n',params.submat(s),submat2(z));
end

end

