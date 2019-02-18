%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bandit: preprocessing pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Requires eeglab

% The pipeline creates intermediate files for each preprocessing step. If
% you change something in a previous step, make sure to rerun ALL
% subsequent steps too to update the files!

%% LOAD DATA AND SET PATHS

clc
clear

% Set path
params.main             = fullfile(''); % change path to folder of analysis scripts
params.toolbox.eeglab   = fullfile(''); % change path to folder of eeglab


if isempty(params.main)
    error('Please change params.main to the location of the folder ''Luyckx_etal_2019''.');
end

if isempty(params.toolbox.eeglab)
    error('Please change params.toolbox.eeglab to the location of your version of eeglab.');
end

params.analysis         = fullfile(params.main,'Analysis');
params.data.main        = fullfile(params.main,'Data');

params.data.behav       = fullfile(params.data.main,'Bandit_behav');
params.data.eeg         = fullfile(params.data.main,'Bandit_EEG_raw');
params.data.saveEEG     = fullfile(params.data.main,'Bandit_EEG_preproc');
params.data.finalEEG    = fullfile(params.data.main,'Bandit_EEG');

params.functions.main  	= fullfile(params.main,'Functions');
params.toolbox.eeglab   = fullfile(params.toolmain, eeglabversion);

cd(params.analysis);
addpath(genpath(params.data.main));
addpath(genpath(params.functions.main));
addpath(params.toolbox.eeglab);

params.whichPhase       = 'test'; % we only analyse the test phase of the bandit task

%% Settings

% Participants to exclude
params.exclude_sub  = [2,8,104];
params.submat       = setdiff([1:24,101:125],params.exclude_sub); % subject indices remaining subjects
params.trialz       = 1:600; % which trials to use?

% Bad channels
switch params.whichPhase
    case 'test'
        
        params.badchan = {
            {'F6','AF3','T8','TP8','TP7'},... % Subject 1
            {'OZ','O2','PO8','AF8','PO7','F6'},... % Subject 3
            {'FT7','T7','TP7','FT8','T8','FP1','AF7'},... % Subject 4
            {},... % Subject 5
            {'FP1','FP2','AF7','F7','AF8','F6','F4'},... % Subject 6
            {'FT7','T8','AF7','F7','FP2','AF8'},... % Subject 7
            {'T7'},... % Subject 9
            {'FP2','O1','PO8','O2','AF7','FP1','AF8'},... % Subject 10
            {'T7','T8'},... % Subject 11
            {'AF7','F4','O2','T8','F1'},... % Subject 12
            {'OZ','O1','PO7','FC4','FC6','AF3','FP1','FP2','AF4','TP7'},... % Subject 13
            {},... % Subject 14
            {'O1','OZ','PO7','PO3','O2','PO4','CP4','P7','CP6','AF8','FT8','F8'},... % Subject 15
            {'T8','C3','CP3','F8'},... % Subject 16
            {'T8','P7','TP7','F2','C2'},... % Subject 17
            {'AF3','C1','P7'},... % Subject 18
            {'O1','PO7','O2','PO8','PO4','T7','FP2'},... % Subject 19
            {'AF7','F5','FP2','AF4','F6','AF8','F4'},... % Subject 20
            {'FC1','FC2','PO4','T8','T7','FT8','TP8','TP7'},... % Subject 21
            {'T7','T8','FP1','FP2','AF8','FC6','TP7','O2','OZ','FC5','C2'},... % Subject 22
            {},... % Subject 23            
            {'OZ','O1','PO7','FT7','T7','T8','O2','FP1'},... % Subject 24
            {'FP1','FP2','FPZ','AF8','F8','C4','FC6','F1','PO8','T7','F6'},... % Subject 101
            {'T7','TP7','O1','FC4','CP3','O2','PO8'},... % Subject 102
            {'T8','TP7','P7'},... % Subject 103
            {'FT8','TP7','F7','FP2','AF3','POZ','P4','CP4'},... % Subject 105
            {'FP1','F6','C1','F5','AF8','FP2','AF7'},... % Subject 106
            {'TP8','AF8','F6','OZ','T7','C4','T8'},... % Subject 107
            {'F6','F4','AF8','AF4','FP2','FP1','AF3','FPZ','T7','CP6'},... % Subject 108
            {'T8','CP2','O2','FC4','P4','C4','PO4','PO8'},... % Subject 109
            {'FT8','P8','F2','P7','TP8','AF3','POZ'},... % Subject 110
            {'C1','FC4','PZ','FP2','AF4','POZ','FT8','CP1','F1','FCZ'},... % Subject 111
            {},... % Subject 112
            {'FP1','AF8','AF7','F6','FP2'},... % Subject 113
            {'P4','P2','FC5','PO3','C2','FC4','CP1','CP3'},... % Subject 114
            {'FP2','AF8','P3','OZ','F4','F6','FC6','AF7','F8','FT8','T8'},... % Subject 115
            {'T8','C4','TP7','CP5','T7','C2','TP8','CP3','F4','P7','P5','FT8'},... % Subject 116
            {'T8','T7','C5','FP2','P8','CP6','F4','F2','FC5','TP7'},... % Subject 117
            {'T8','FCZ','FP2','FP1','FC4'},... % Subject 118
            {'O2','O1','PO8','FP2','PO7','AF8'},... % Subject 119
            {'FC2','PO3','CP1','FT8'},... % Subject 120
            {'CP1','PO3','O1','T7','CPZ','PZ'},... % Subject 121
            {'CP1','C2'},... % Subject 122
            {'T8','FC1','P4','P6','C1','C3','FP1','FP2','FPZ'},... % Subject 123
            {'FP1','AF8'},... % Subject 124
            {'OZ','PZ','AF8','F8','O2'},... % Subject 125
            };
end

params.resamp           = 250; % resampling rate
params.filter.highpass 	= 1; % high-pass filter
params.filter.lowpass 	= 40; % low-pass filter
params.triggerz.all     = {1 2 11 12 21 22 31 33 41 42};

params.triggerz.long    = params.triggerz.all(1);
params.triggerz.short   = params.triggerz.all([3,5]);

params.epoch.long       = [-.5 3]; % time window of long epoch (check what your first trigger will be!)
params.epoch.short      = [-.25 .75]; % time window of short epoch
params.baseline.long  	= [-500 0]; % baseline trial epoch
params.baseline.short   = [-250 0]; % baseline sample epoch / make empty if you don't want it

%% Steps

do.parpooling       = false;

do.import           = false;    % create EEGLAB files/downsample/remove mastoid A1, EOG
do.filtering        = false;    % filter data
do.badchan          = false;     % remove bad channels and interpolate
do.averageref       = false;     % average referencing
do.longepoching     = false;     % extract long epochs and baseline correct
do.artefact         = false;     % then we need to visually inspect the data and reject epochs for each subject
do.rejectagain      = false;    % when you changed stuff in the previous steps, but don't want to reject trials manually again
do.ica              = false;    % use automatic ICA (can be quite long)
do.componentrej     = false;    % component rejection of ICA
do.shortepoching    = false;    % epoch on sample level and possible re-baseline
do.extract          = false;    % make mat-file for data analysis
do.extractShort     = false;    % extract short (true) or long (false) epoch?

%% Supercomputer settings

if do.parpooling
    numWorkers = length(params.submat);
    parpool(length(params.submat));
else
    numWorkers = 0;
end

%% Preprocessing pipeline
try
    
    %% 1. CREATE EEGLAB FILES and DOWNSAMPLE DATA
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        inputfile   = sprintf('Donkey_%03d_%s.dat',params.submat(s),params.whichPhase);
        outputfile  = sprintf('Donkey_%03d_%s_downsampled.set',params.submat(s),params.whichPhase);
        
        fEEG_import(params,do,s,inputfile,outputfile);
    end
    
    %% 2. FILTER DATA
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        inputfile   = sprintf('Donkey_%03d_%s_downsampled.set',params.submat(s),params.whichPhase);
        outputfile 	= sprintf('Donkey_%03d_%s_filter.set',params.submat(s),params.whichPhase);
        
        fEEG_filter(params,do,s,inputfile,outputfile);
    end
    
    %% 3. BAD CHANNEL DETECTION AND INTERPOLATION
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        inputfile 	= sprintf('Donkey_%03d_%s_filter.set',params.submat(s),params.whichPhase);
        outputfile 	= sprintf('Donkey_%03d_%s_interpolated.set',params.submat(s),params.whichPhase);
        
        fEEG_badchan(params,do,s,inputfile,outputfile);
    end
    
    %% 4. AVERAGE REFERENCING
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        inputfile   = sprintf('Donkey_%03d_%s_interpolated.set',params.submat(s),params.whichPhase);
        outputfile 	= sprintf('Donkey_%03d_%s_avref.set',params.submat(s),params.whichPhase);
        
        fEEG_averageref(params,do,s,inputfile,outputfile);
    end
    
    %% 5. LONG EPOCHING
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        inputfile   = sprintf('Donkey_%03d_%s_avref.set',params.submat(s),params.whichPhase);
        outputfile 	= sprintf('Donkey_%03d_%s_longepoch.set',params.submat(s),params.whichPhase);
        
        fEEG_longepoch(params,do,s,inputfile,outputfile);
    end
    
    %% 6. REJECT ARTEFACTS (partially manually)
    
    if do.artefact
        
        % click on bad trials + UPDATE MARKS before closing plot window!!
        % Press any key in command window to continue automatic rejection
        % of marked epochs
        
        subject = [];
        while isempty(subject)
            subject = input('Which subject? ');
            if ~any(subject == params.submat)
                error('Subject index not found.');
            end
        end
        
        inputfile   = sprintf('Donkey_%03d_%s_longepoch.set',subject,params.whichPhase);
        outputfile	= sprintf('Donkey_%03d_%s_rejtrials.set',subject,params.whichPhase);
        mfile       = sprintf('Donkey_%03d_%s_rejectedTrials.mat',subject,params.whichPhase);
        
        if ~isempty(subject)
            disp(' ');
            disp(['Rejecting trials of subject ' num2str(subject)])
            
            % restart eeglab (because of memory issues)
            close all
            clear ALLEEG EEG CURRENTSET ALLCOM
            [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
            
            % load dataset
            EEG = pop_loadset(inputfile, params.data.saveEEG);
            [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, 1);
            
            % Open EEG plot
            pop_eegplot(EEG,1,1,0);
            pause
            EEG = ALLEEG(2);
            
            % Reject marked epochs
            marksUpdated = questdlg('Done with updating marks?');
            if strmatch(marksUpdated, 'Yes')
                
                % Save rejected trials in m-file
                clear rejectedTrialz
                rejectedTrialz = EEG.reject.rejmanual;
                save(fullfile(params.data.finalEEG,mfile),'rejectedTrialz');
                fprintf('\nRejected trials m-file subject %d saved.\n',subject);
                
                % Reject artefacts
                EEG = pop_rejepoch(EEG,rejectedTrialz);
                
                % Save file
                EEG.setname     = outputfile;
                EEG             = pop_saveset(EEG, outputfile, params.data.saveEEG);
                
                fprintf('\nRejected trials subject %d finished and saved.\n',subject);
                
            end
        end
    end
    
    %% EXTRA: reject trials again
    % After you've done something in the previous steps, but you don't want to
    % manually start rejecting all the trials you've rejected before.
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        inputfile   = sprintf('Donkey_%03d_%s_longepoch.set',params.submat(s),params.whichPhase);
        outputfile	= sprintf('Donkey_%03d_%s_rejtrials.set',params.submat(s),params.whichPhase);
        mfile       = sprintf('Donkey_%03d_%s_rejectedTrials.mat',params.submat(s),params.whichPhase);
        
        fEEG_rejectagain(params,do,s,inputfile,outputfile,mfile);
    end
    
    %% 7. ICA
    
    if do.ica
        cstart = clock;
    end
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        inputfile   = sprintf('Donkey_%03d_%s_rejtrials.set',params.submat(s),params.whichPhase);
        outputfile 	= sprintf('Donkey_%03d_%s_ica.set',params.submat(s),params.whichPhase);
        
        fEEG_ica(params,do,s,inputfile,outputfile);
    end
    
    if do.ica
        cfinish = clock;
        fprintf('\nICA started at %d:%d and finished at %d:%d.\n',cstart(4),cstart(5),cfinish(4),cfinish(5));
    end
    
    %% 8. REMOVE ICA COMPONENTS
    
    if do.componentrej
        
        close all
        clear ALLEEG EEG CURRENTSET ALLCOM
        
        subject = [];
        while isempty(subject)
            subject = input('Which subject? ');
            if ~any(subject == params.submat)
                error('Subject index not found.');
            end
        end
        
        inputfile   = sprintf('Donkey_%03d_%s_ica.set',subject,params.whichPhase);
        outputfile	= sprintf('Donkey_%03d_%s_prunedica.set',subject,params.whichPhase);
        mfile       = sprintf('Donkey_%03d_%s_rejectedComps.mat',subject,params.whichPhase);
        
        if ~isempty(subject)
            disp(' ');
            disp(['Removing components of subject ' num2str(subject)])
            
            % restart eeglab (because of memory issues)
            [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
            
            % load dataset
            EEG = pop_loadset(inputfile, params.data.saveEEG);
            [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, 1);
            
            % Open EEG plot
            pop_selectcomps(EEG,1:35); % plot the component map for selection
            pop_eegplot(EEG, 0, 1, 0); % plot ICA compoments timecourse % change scale to 10 and set full screen
            pause
            EEG = ALLEEG(2);
            
            % Reject marked epochs
            compsUpdated = questdlg('Done with marking components?');
            if strmatch(compsUpdated, 'Yes')
                
                % Save removed components in m-file
                clear rejectedCompz
                rejectedCompz = find(EEG.reject.gcompreject == 1);
                save(fullfile(params.data.saveEEG,mfile),'rejectedCompz');
                disp(' ');
                disp(['Removed components m-file subject ' num2str(subject) ' saved.']);
                fprintf('\nRemoved components m-file subject %d saved.\n',subject);
                
                % Remove components
                EEG = pop_subcomp(EEG);
                
                % Save file
                EEG.setname     = outputfile;
                EEG             = pop_saveset(EEG, outputfile, params.data.saveEEG);
                
                fprintf('\nRemoving %d components from subject %d finished and saved.\n',length(rejectedCompz),subject);
                
            end
        end
    end
    
    %% 9. SHORT EPOCH
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        inputfile   = sprintf('Donkey_%03d_%s_prunedica.set',params.submat(s),params.whichPhase);
        outputfile	= sprintf('Donkey_%03d_%s_shortepoched.set',params.submat(s),params.whichPhase);
        
        fEEG_shortepoch(params,do,s,inputfile,outputfile);
    end
    
    %% 12. EXTRACT RELEVANT DATA IN MAT-FILE
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        if do.extractShort
            inputfile   = sprintf('Donkey_%03d_%s_shortepoched.set',params.submat(s),params.whichPhase);
            outputfile	= sprintf('Donkey_%03d_%s_EEG_samples.mat',params.submat(s),params.whichPhase);
        else
            inputfile   = sprintf('Donkey_%03d_%s_prunedica.set',params.submat(s),params.whichPhase);
            outputfile	= sprintf('Donkey_%03d_%s_EEG_trials.mat',params.submat(s),params.whichPhase);
        end
        
        fEEG_extractmat(params,do,s,inputfile,outputfile);
        
    end
    
    %% End parpool session
    
    if do.parpooling
        delete(gcp());
    end
    
catch ME
    
    if do.parpooling
        delete(gcp());
    end
    
    rethrow(ME)
    disp(' ');
    disp('Try loop failed.');
    
end