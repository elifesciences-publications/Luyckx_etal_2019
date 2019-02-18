%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Donkey_main: runner
%
% All functions written by Fabrice Luyckx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

execute = questdlg('Clear all and start experiment?');
%execute = 'Yes';

if strmatch(execute, 'Yes')
    
    clc
    clear 
    close all
    
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
    datafolder = fullfile(datapath,'/Donkey_data/');
    
    addpath(genpath(fullfile(datapath,'Donkey_functions/'))); % add function folder with subfolder
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
    nsamp       = 6;                            % number of bandits
    aborted     = 0;
    
    randomise   = 1;                            
    if randomise == 1 % variable to activate randomisation
        % Seed random number generator
        rng('shuffle');
    end   
    
    %% Image mapping
    if randomise == 1
        imgMap = Shuffle(1:nsamp);
    else
        imgMap = 1:nsamp;        
    end
    
    %% Click task
    whichPhase      = 1;
    clickTrialz     = 36*nsamp;
    clickBlockz     = 1;    

    disp('Click task initiated.')
    [data, stim, time, triggers, aborted] = Donkey_click(inEEG,randomise,imgMap,ppnr,clickTrialz,clickBlockz,nsamp);
    
    if aborted == 0
        
        % Save data
        participant.filename = sprintf('Donkey_ppt_%s_%3s_click.mat',participant.name,datestr(now,'yyyymmddHHMMSS'));
        datafile = fullfile(datapath,'/Donkey_data/',participant.filename);
        save(datafile,'participant','data','stim','time','triggers');
        
        Donkey_send_email;
        clear data stim time
        
    end
    
    %% Learning phase    
    if aborted == 0
        whichPhase      = 2;
        practTrialz     = 120;
        practBlockz     = 2;
        
        disp('Learning phase initiated.')
        [data, stim, time, triggers, aborted] = Donkey_trial(inEEG,randomise,whichPhase,imgMap,ppnr,practTrialz,practBlockz,nsamp);
        
        % Save data and send e-mail
        if aborted == 0
           
            % Save data
            participant.filename = sprintf('Donkey_ppt_%s_%3s_learn.mat',participant.name,datestr(now,'yyyymmddHHMMSS'));
            datafile = fullfile(datapath,'/Donkey_data/',participant.filename);
            save(datafile,'participant','data','stim','time','triggers');
            
            % Send email
            Donkey_send_email;
            clear data stim time
        end        
    end
    
    %% Test phase
    if aborted == 0
        whichPhase      = 3;
        testTrialz      = 600;                          % number of trials
        testBlockz      = 10;                            % number of blocks
        
        disp('Testing phase initiated.')
        [data, stim, time, triggers, aborted] = Donkey_trial(inEEG,randomise, whichPhase,imgMap, ppnr, testTrialz, testBlockz,nsamp);
        
        % Send e-mail with data
        if aborted == 0
            % Save data
            participant.filename = sprintf('Donkey_ppt_%s_%3s_test.mat',participant.name,datestr(now,'yyyymmddHHMMSS'));
            datafile = fullfile(datapath,'/Donkey_data/',participant.filename);
            save(datafile,'participant','data','stim','time','triggers');
            
            Donkey_send_email;
        end
    end
    
end

%% EXPLANATION DIFFERENT FUNCTIONS
%
% Donkey_initialise = define all variables for data and stim(ulus) structures,
% very handy to check when something went wrong in the randomisation or
% trial files.
%
% Donkey_setup = give all variables the real values and option to randomise
% everything. Use this to check whether randomisation does what it's
% supposed to do and check underlying distributions.
%
% Donkey_trial = run the real experiment without all the fuzz around it
%
% Donkey_main = file that executes the whole experiment. Here you define the
% actual values of the number of trials etc. in your experiment.













