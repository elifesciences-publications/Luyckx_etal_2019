%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Numbers: EEG preprocessing runner
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

params.data.behav       = fullfile(params.data.main,'Numbers_behav');
params.data.eeg         = fullfile(params.data.main,'Numbers_EEG_raw');
params.data.saveEEG     = fullfile(params.data.main,'Numbers_EEG_preproc');
params.data.finalEEG    = fullfile(params.data.main,'Numbers_EEG');

params.functions.main  	= fullfile(params.main,'Functions');

cd(params.analysis);
addpath(genpath(params.data.main));
addpath(genpath(params.functions.main));
addpath(params.toolbox.eeglab);

%% Settings

% Participants to exclude
params.exclude_sub      = [2,8,104];
params.submat           = setdiff([1:24,101:125],params.exclude_sub); % subject indices

params.trialz           = 1:300; % which trials to use?
params.resamp           = 250; % resampling rate
params.filter.highpass 	= 1; % high-pass filter
params.filter.lowpass 	= 40; % low-pass filter
params.triggerz.all     = {1 10 11 12 13 14 15 16 17 18 19 31 33 41 42};

for i = 1:50
    params.triggerz.long{i} = 100+(i-1);
end
params.triggerz.short   = params.triggerz.all(2:11);

params.epoch.long       = [-1 5.5]; % time window of long epoch (check what your first trigger will be!)
params.epoch.short      = [-.065 .85]; % time window of short epoch
params.baseline.long  	= [-1000 0]; % baseline trial epoch
params.baseline.short   = [-65 0]; % baseline sample epoch / make empty if you don't want it

% Bad channels
params.badchan = {
    {'F6','AF3','T8','TP8'},... % Subject 1
    {'OZ','O2','PO8','AF8','CP5','P5','AF7'},... % Subject 3
    {'FT7','T7','TP7','TP8','P8','FT8','T8','FP1'},... % Subject 4
    {'T7'},... % Subject 5
    {'FP1','FP2','AF7','AF8','O2','PO8','F7','F2'},... % Subject 6
    {'AF7','FT8'},... % Subject 7
    {'T7','FP1','AF7'},... % Subject 9
    {'FP2','FP1','FPZ','AF8','F5','AF7'},... % Subject 10
    {'T7','T8','FT8'},... % Subject 11
    {'AF7','FP2','FP1','AF8','OZ','O1','T8'},... % Subject 12
    {'O1','PO7','AF8','AF3','OZ','TP7'},... % Subject 13
    {},... % Subject 14
    {'O1','O2','OZ','CP6','PO3','PO4','PO8','P7','PO7','C4','FCZ','P6'},... % Subject 15
    {'PO4','P8','PO8','O2','FP2'},... % Subject 16
    {'TP7','P7','T8','C1','O1','OZ','PO7'},... % Subject 17
    {'T8','FP2'},... % Subject 18
    {'PO7','O1','O2','PO8','P5'},... % Subject 19
    {'AF7','F5','FP2','AF8','F6','F4','P3'},... % Subject 20
    {'T8','T7','PO4','PO8','OZ','O1'},... % Subject 21
    {'T7','T8','FC6','O2'},... % Subject 22
    {'T7','T8','TP7','P4'},... % Subject 23
    {'OZ','O1','PO7','T7','T8','FP1','P7'},... % Subject 24
    {'FPZ','FP1','FP2','T7','AF8','F8','F6','FC6','AF3','TP8'},... % Subject 101
    {'T7','PO3','O1','FP2','PO8','FC2','F2','AF7','FC5','FC3'},... % Subject 102
    {'T8','FP2','AF4','F4','F6','P8','TP7'},... % Subject 103
    {'T7','TP7','AF8','F7'},... % Subject 105
    {'FP1','FP2','AF8','F6','F5','AF4'},... % Subject 106
    {'T7','C2','O1','F6','FP1','FP2','AF8','FC5','TP8','T8','CPZ'},... % Subject 107
    {'FPZ','FP1','AF3','T7','FP2','AF4','AF8','F6','T8','CP1','F8','F4'},... % Subject 108
    {'T8','C6','C4','P6','AF4','FC4','P4','FC6'},... % Subject 109
    {'F2'},... % Subject 110
    {'F5','F3','FC5','FC3','O2','CP4','FT8','F1','C5','PZ','FC1'},... % Subject 111
    {},... % Subject 112
    {'FP1','FP2','AF7','AF8','F6','F2','C2'},... % Subject 113
    {'AF8','P2','F6','P4','F2','P8'},... % Subject 114
    {'T8','FT8','FP2','OZ','AF8','AF4','F8','AF7','F4','F6','FC4','FC6'},... % Subject 115
    {'CP5','F3','FPZ','T7','OZ'},... % Subject 116
    {'T7','TP7','T8','FC5','FP1','FP2'},... % Subject 117
    {'T8','AF3','FP1','FP2'},... % Subject 118
    {'C1','C4','O1','O2','FP2','FP1','PO8','PO7'},... % Subject 19
    {'O1','PO4','T7','F6','FT7'},... % Subject 120
    {'P1','O2','AF7','P3','CP2'},... % Subject 121
    {'TP8','P4','CP1','CZ'},... % Subject 122
    {'T8','FP1','AF3','FT8','AF8','C6','FCZ','FPZ','FP2'},... % Subject 123
    {'AF8','FP1','FT8'},... % Subject 124
    {'AF8','AF7','PO3','F6','FC2','F8','FP2','F4','P5','TP7'},... % Subject 125
    };

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
        
        inputfile   = sprintf('Numbers_%03d.dat',params.submat(s));
        outputfile  = sprintf('Numbers_%03d_downsampled.set',params.submat(s));
        
        fEEG_import(params,do,s,inputfile,outputfile);
    end
    
    %% 2. FILTER DATA
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        inputfile   = sprintf('Numbers_%03d_downsampled.set',params.submat(s));
        outputfile 	= sprintf('Numbers_%03d_filter.set',params.submat(s));
        
        fEEG_filter(params,do,s,inputfile,outputfile);
    end
    
    %% 3. BAD CHANNEL DETECTION AND INTERPOLATION
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        inputfile 	= sprintf('Numbers_%03d_filter.set',params.submat(s));
        outputfile 	= sprintf('Numbers_%03d_interpolated.set',params.submat(s));
        
        fEEG_badchan(params,do,s,inputfile,outputfile);
    end
    
    %% 4. AVERAGE REFERENCING
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        inputfile   = sprintf('Numbers_%03d_interpolated.set',params.submat(s));
        outputfile 	= sprintf('Numbers_%03d_avref.set',params.submat(s));
        
        fEEG_averageref(params,do,s,inputfile,outputfile);
    end
    
    %% 5. LONG EPOCHING
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        inputfile   = sprintf('Numbers_%03d_avref.set',params.submat(s));
        outputfile 	= sprintf('Numbers_%03d_longepoch.set',params.submat(s));
        
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
        
        inputfile   = sprintf('Numbers_%03d_longepoch.set',subject);
        outputfile	= sprintf('Numbers_%03d_rejtrials.set',subject);
        mfile       = sprintf('Numbers_%03d_rejectedTrials.mat',subject);
        
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
                rejectedTrialz = EEG.reject.rejmanual(params.trialz);
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
    
    for s = 1:length(params.submat)
        
        inputfile   = sprintf('Numbers_%03d_longepoch.set',params.submat(s));
        outputfile	= sprintf('Numbers_%03d_rejtrials.set',params.submat(s));
        mfile       = sprintf('Numbers_%03d_rejectedTrials.mat',params.submat(s));
        
        fEEG_rejectagain(params,do,s,inputfile,outputfile,mfile);
    end
    
    %% 7. ICA
    
    if do.ica
        cstart = clock;
    end
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        inputfile   = sprintf('Numbers_%03d_rejtrials.set',params.submat(s));
        outputfile 	= sprintf('Numbers_%03d_ica.set',params.submat(s));
        
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
        
        inputfile   = sprintf('Numbers_%03d_ica.set',subject);
        outputfile	= sprintf('Numbers_%03d_prunedica.set',subject);
        mfile       = sprintf('Numbers_%03d_rejectedComps.mat',subject);
        
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
        
        inputfile   = sprintf('Numbers_%03d_prunedica.set',params.submat(s));
        outputfile	= sprintf('Numbers_%03d_shortepoched.set',params.submat(s));
        
        fEEG_shortepoch(params,do,s,inputfile,outputfile);
    end
    
    %% 12. EXTRACT RELEVANT DATA IN MAT-FILE
    
    parfor (s = 1:length(params.submat),numWorkers)
        
        if do.extractShort
            inputfile   = sprintf('Numbers_%03d_shortepoched.set',params.submat(s));
            outputfile	= sprintf('Numbers_%03d_EEG_samples.mat',params.submat(s));
        else
            inputfile   = sprintf('Numbers_%03d_prunedica.set',params.submat(s));
            outputfile	= sprintf('Numbers_%03d_EEG_trials.mat',params.submat(s));
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