%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Numbers_main: runner
%
% All functions written by Fabrice Luyckx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

execute = questdlg('Clear all and start experiment?');
%execute = 'Yes';

if strmatch(execute, 'Yes')
    
    clc
    clear 
    close all
    
    % Using EEG?
    EEGquest   = questdlg('Using EEG system?');
    switch EEGquest
        case 'Yes'
            inEEG = true;
        case 'No'
            inEEG = false;
        case 'Cancel'
            return;
    end
    
    %% Set path
    datapath = cd;
    datafolder = fullfile(datapath,'/Numb_data/');
    
    addpath(genpath(fullfile(datapath,'Numb_functions/'))); % add function folder with subfolder
    addpath(datafolder); % add datafolder
    
%% Create participant structure

    argindlg = inputdlg({'Participant number   ','Gender (M/F)','Age','Hand (L/R)','Which part?'},'',1,{'000','','','R',''});
    if isempty(argindlg)
        return;
    else
        participant                 = struct;
        participant.name            = upper(argindlg{1});
        participant.gender          = argindlg{2};
        participant.age             = argindlg{3};
        participant.handedness      = argindlg{4};
        participant.part            = argindlg{5};
    end
    
    %% Get OS
    os = computer;
    
    %% Set keys
    KbName('UnifyKeyNames');
    
    %% What I need to know
    ppnr        = str2num(participant.name);    % participant number
    ntrials     = 300;                          % number of trials % should be 300
    nblocks     = 6;                            % number of blocks
    nsamp       = 10;                            % number of samples
    numrange    = 6;                            % range of numbers
    
    randomise   = 1;                            
    if randomise == 1 % variable to activate randomisation
        % Seed random number generator
        rng('shuffle');
    end   
    
    %% Demo
    practice = 2;
    demoTrialz = 4;
    demoBlockz = 1;
    
    disp ('Demo initiated.')
    [~, ~, ~, ~, aborted] = Exp_trial(inEEG,randomise, practice, ppnr, demoTrialz,demoBlockz,nsamp,numrange);
    
    %% Practice trials
    if aborted == 0
        practice = 1;
        practTrialz = 10;
        practBlockz = 1;
        
        disp('Practice initiated.')
        [~, ~, ~, ~, aborted] = Exp_trial(inEEG,randomise,practice,ppnr,practTrialz,practBlockz,nsamp,numrange);
    end   
    
    %% Run experiment
    if aborted == 0
        practice = 0;
        disp('Experiment initiated.')
        [data, stim, time, triggers, aborted] = Exp_trial(inEEG,randomise, practice, ppnr, ntrials, nblocks, nsamp, numrange);
    end
    
    if aborted == 0
        %% Save data
        participant.filename = sprintf('Numb_ppt_%s_%3s.mat',participant.name,datestr(now,'yyyymmddHHMMSS'));
        datafile = fullfile(datapath,'/Numb_data/',participant.filename);
        save(datafile,'participant','data','stim','time','triggers');
        
        %% Send e-mail with data
        Exp_send_email;
    end
    
end

%% EXPLANATION DIFFERENT FUNCTIONS
%
% Exp_initialise = define all variables for data and stim(ulus) structures,
% very handy to check when something went wrong in the randomisation or
% trial files.
%
% Exp_setup = give all variables the real values and option to randomise
% everything. Use this to check whether randomisation does what it's
% supposed to do and check underlying distributions.
%
% Exp_trial = run the real experiment without all the fuzz around it
%
% Exp_main = file that executes the whole experiment. Here you define the
% actual values of the number of trials etc. in your experiment.













