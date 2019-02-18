%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 1C: RDM individual tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For Figure 1 - supplement 1B, change params.disttype to 'pearson'.

%% LOAD DATA

clc
clear

% Paths
savefolder          = '';
figfolder           = 'RSA_indiv'; % folder to save figures to
params.whichPhase   = 'test'; % use test phase data
params.disttype     = 'euclidean'; % Fig1C -> 'euclidean', Fig1-suppl1B -> 'pearson'

% Load stuff
Config_plot; % load plot variables

% Logicals
do.RDM              = true; % create RDM data, if it doesn't already exist
do.saveRDM          = true; % also save the RDM data, if it doesn't already exist

do.smooth           = true; % smooth the data?
do.plotting         = true; % plot?
do.signif           = true; % test significance?

%% Extra variables

% Load RL fits
Bandit_load;
load(fullfile(paths.data.model,'Modelfit_full_test_RL'));

% Replace bandit index with perceived index
subjranks = 0.*stim.combo;

for t = 1:params.ttrials
    for i = 1:params.nsamp
        subjranks(t,i) = find(stim.combo(t,i) == mod.ranks(t,:));
    end
end

% Task names
params.tasknames    = {'numbers','bandit'};
condnamez           = {'1','2','3','4','5','6';'b1','b2','b3','b4','b5','b6'};

%% Get mutual time

mutual_time     = -64:4:748;
ntimepoints     = length(mutual_time);

%% Concatenate EEG RDM

for c = 1:2 % loop through tasks
    for s = 1:params.nsubj % subject correlations
        
        if c == 1
            fprintf('\nNumerical: correlating model - EEG RDM subject %d.\n',params.submat(s));
            
            % Load paths/data
            if s == 1
                Numbers_load;
            end
            
            % File to load
            inputfile   = sprintf('Numbers_%03d_EEG_RDM_%s.mat',params.submat(s),params.disttype);
            
        elseif c == 2
            fprintf('\nBandit: correlating model - EEG RDM subject %d, %s phase.\n',params.submat(s),params.whichPhase);
            
            % Load paths/data
            if s == 1
                Bandit_load;
            end
            
            % File to load
            inputfile  = sprintf('Donkey_%03d_%s_EEG_RDM_%s',params.submat(s),params.whichPhase,params.disttype);
        end
        
        % Load RDM
        load(fullfile(paths.data.EEG.RDM,inputfile));
        
        % Smooth data
        if do.smooth
            wdwsz    = 60/4; % size convolution kernel (ms/downsampling rate)
            rdm.data = smoothRDM(rdm.data,wdwsz);
        end
        
        % Initialise and store time points (different sizes of epochs)
        if s == 1
            if c == 1
                eegrdm = zeros(6,6,ntimepoints,params.nsubj,2);
                timepoints1 = rdm.timepoints;
                eeg1wdw 	= timepoints1 >= mutual_time(1) & timepoints1 <= mutual_time(end); % numbers idx
            elseif c == 2
                timepoints2 = rdm.timepoints;
                eeg2wdw 	= timepoints2 >= mutual_time(1) & timepoints2 <= mutual_time(end); % bandit idx
            end
        end
        
        % Concatenate EEG RDM
        if c == 1
            eegrdm(:,:,:,s,c) = rdm.data(:,:,eeg1wdw);
        elseif c == 2
            eegrdm(:,:,:,s,c) = rdm.data(:,:,eeg2wdw);
        end
        
    end

end

%% EEG RDM

for c = 1:2

    timewdw = [200 700];
    timeidx = mutual_time >= timewdw(1) & mutual_time <= timewdw(end);
    plotrdm = squeeze(mean(mean(eegrdm(:,:,timeidx,:,c),3),4));
    
    maplims = [min(plotrdm(plotrdm ~= 0)) max(max(plotrdm(:)))];
    
    figure;
    imagesc(plotrdm,maplims);
    colormap(hot);
    axis square
    ax  = gca;
    cax = colorbar('southoutside');
    cax.Limits = round(maplims,2);
    set(ax,'XTick',1:length(condnamez(c,:)),'XTickLabel',condnamez(c,:));
    set(ax,'YTick',1:length(condnamez(c,:)),'YTickLabel',condnamez(c,:));
    set(cax,'YTick',round(maplims,2));
    ax.XAxis.FontSize   = 24;
    ax.YAxis.FontSize   = 24;
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontWeight = 'bold';
    cax.FontSize = 24;
    
end