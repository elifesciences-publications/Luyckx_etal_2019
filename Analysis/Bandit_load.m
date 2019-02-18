%% Load data/set paths for bandit task

% Adjust these paths
paths.main                  = fullfile(''); % change path to location of folder
paths.toolbox.eeglab        = fullfile(''); % path to location of eeglab version

if isempty(paths.main)
    error('Please change paths.main to the location of the folder ''Luyckx_etal_2019''.');
end

if isempty(paths.toolbox.eeglab)
    error('Please change paths.toolbox.eeglab to the location of your version of eeglab.');
end

% Set path
paths.analysis          = fullfile(paths.main,'Analysis');

paths.data.main         = fullfile(paths.main,'Data','Bandit');
paths.data.behav        = fullfile(paths.data.main,'Bandit_behav');
paths.data.model        = fullfile(paths.data.main,'Bandit_model');
paths.data.EEG.raw      = fullfile(paths.data.main,'Bandit_EEG');
paths.data.EEG.RDM      = fullfile(paths.data.main,'Bandit_RDM');
paths.data.save         = fullfile(paths.data.main,savefolder);

paths.functions.main    = fullfile(paths.main,'Functions');
paths.figures.current   = fullfile(paths.main,'Figures',figfolder);

% Create folder for figures if doesn't exist yet
if exist(paths.figures.current,'dir') ~= 7
    mkdir(fullfile(paths.main,'Figures'),figfolder);
end

cd(paths.analysis);
addpath(genpath(paths.data.main));
addpath(genpath(paths.functions.main));
addpath(paths.toolbox.eeglab);

% Load behavioural data
load(fullfile(paths.data.behav,sprintf('Donkey_fulldata_%s_all',params.whichPhase)));

% Load channel info
load(fullfile(paths.data.EEG.raw,'chanlocs_curry.mat'));

%% Variables

participant.part    = str2num(participant.part);

params.submat       = unique(data.sub)';
params.nsubj        = length(params.submat);
params.ttrials      = length(data.sub);
params.ntrials      = length(data.sub(data.sub == params.submat(1)));
params.nblocks      = length(unique(data.block));
params.btrials      = params.ntrials/params.nblocks;
params.nstim        = size(stim.proba,2);
params.nsamp        = size(stim.combo,2);
params.nparts       = length(unique(data.part));

for z = 1:params.nsubj
    params.frame(z) = unique(data.frame(data.sub == params.submat(z)));
end
