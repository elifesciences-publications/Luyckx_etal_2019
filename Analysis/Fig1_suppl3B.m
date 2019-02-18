%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 1 - suppl3B: RDM first vs second bandit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DATA

clc
clear

% Paths
savefolder          = 'Bandit_RDM';
figfolder           = 'RSA_indiv'; % folder to save figures to
params.whichPhase   = 'test'; % use test phase data
params.disttype     = 'euclidean';

% Load stuff
Bandit_load;
Config_plot; % load plot variables

% Logicals
do.RDM              = false; % create RDM data, if it doesn't already exist
do.saveRDM          = false; % also save the RDM data, if it doesn't already exist

do.smooth           = true; % smooth the data?
do.signif           = true; % test significance?

%% Extra variables

% Load RL fits
load(fullfile(paths.data.model,'Modelfit_full_test_RL'));

% Replace bandit index with perceived index
subjranks = 0.*stim.combo;

for t = 1:params.ttrials
    for i = 1:params.nsamp
        subjranks(t,i) = find(stim.combo(t,i) == mod.ranks(t,:));
    end
end

%% Create RDM

if do.RDM
    for s = 1:params.nsubj
        
        inputfile  	= sprintf('Donkey_%03d_%s_',params.submat(s),params.whichPhase);
        outputfile  = sprintf('Donkey_%03d_%s_EEG_RDM_%s',params.submat(s),params.whichPhase,params.disttype);
        
        data.conds  = subjranks;
        
        CreateRDM_bandits(s,inputfile,outputfile,do,data,params,paths.data.EEG);
    end
end


%% Plot grand average neural RDMs

for b = 1:2
    
    for s = 1:params.nsubj
        
        fprintf('\nConcatenate subject %d.\n',params.submat(s));
        
        % Load data
        inputfile = sprintf('Donkey_%03d_%s_EEG_RDM_%s_bandit%d',params.submat(s),params.whichPhase,params.disttype,b);
        load(fullfile(paths.data.save,inputfile));
        
        if s == 1
            eegrdm = zeros(6,6,size(rdm.data,3),params.nsubj);
        end
        
        % Smoothing RDM
        if do.smooth
            wdwsz  	= 60/4;
            rdm.data = smoothRDM(rdm.data,wdwsz);
        end
        
        eegrdm(:,:,:,s) = rdm.data;
        
    end
    
    %% Plot RDMs
        
    figA    = figure;
    
    timewdw     = [200 700];
    timeidx     = rdm.timepoints >= timewdw(1) & rdm.timepoints <= timewdw(end);
    condnamez   = {'b1','b2','b3','b4','b5','b6'};
    
    plotrdm = squeeze(mean(mean(eegrdm(:,:,timeidx,:),3),4));
    maplims = [min(plotrdm(plotrdm ~= 0)) max(max(plotrdm(:)))];
    
    imagesc(plotrdm,maplims); hold on;
    colormap(hot);
    axis square
    ax  = gca;
    cax = colorbar('southoutside');
    cax.Limits = round(maplims,2);
    set(ax,'XTick',1:length(condnamez),'XTickLabel',condnamez);
    set(ax,'YTick',1:length(condnamez),'YTickLabel',condnamez);
    set(cax,'YTick',round(maplims,2));
    ax.XAxis.FontSize   = 20;
    ax.YAxis.FontSize   = 20;
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontWeight = 'bold';
    cax.FontSize = 20;
        
end







