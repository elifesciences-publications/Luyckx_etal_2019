%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 1 - suppl3A: RSA magnitude first vs second bandit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
do.RDM              = true; % create RDM data, if it doesn't already exist
do.saveRDM          = true; % also save the RDM data, if it doesn't already exist

do.modelcorr        = true; % correlate RDM with model
do.smooth           = true; % smooth the data?
do.plotting         = true; % plot?
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

%% Cross-validation

if do.modelcorr
    
    clear modcont models
    
    % Obtain model RDMs
    mods    = ModelRDM;
    models  = fieldnames(mods);
    
    for b = 1:2 % per bandit
        for s = 1:params.nsubj % Subject correlations
            
            fprintf('\nCorrelating model - EEG RDM subject %d, %s phase.\n',params.submat(s),params.whichPhase);
            
            % Load data
            inputfile = sprintf('Donkey_%03d_%s_EEG_RDM_%s_bandit%d',params.submat(s),params.whichPhase,params.disttype,b);
            load(fullfile(paths.data.save,inputfile));
            
            if s == 1 && b == 1
                rsa_results = zeros(params.nsubj,size(rdm.data,3),2);
            end
            
            % Smoothing RDM
            if do.smooth
                wdwsz    = 60/4; % size convolution kernel
                rdm.data = smoothRDM(rdm.data,wdwsz);
            end
            
            % EEG - Model correlation
            actmod = squareform(mods.(models{1}));
            for t = 1:length(rdm.timepoints)
                rsa_results(s,t,b) = rankCorr_Kendall_taua(actmod,squareform(rdm.data(:,:,t)));
            end
        end
    end
end

%% Plot time series

if do.plotting
    
    mutual_time = [-64 748];
    eegwdw      = rdm.timepoints >= mutual_time(1) & rdm.timepoints <= mutual_time(end); % bandit idx
    tp          = rdm.timepoints(eegwdw);
    plotdat     = rsa_results(:,eegwdw,:);
    
    % Settings
    nmods       = size(plotdat,3);
    modlabels   = {'First bandit','Second bandit'};
    banditcolz  = [colz(2,:).*.8;colz(2,:).*.5];
    
    ylims       = [-.1 .3];
    sigpos      = [-.08, -.07]; % Y position of significance
    
    minX        = min(mutual_time);
    maxX        = max(mutual_time);
    minY        = ylims(1);
    maxY        = ylims(2);
    
    % Figure for line plot
    figR = figure; hold on;
    
    % Plot lines
    plot([minX maxX],[0 0],'k--','LineWidth',1.5);
    plot([0 0],[minY*1.1 maxY*1.1],'k-','LineWidth',1.5);
    steh = myfillsteplot(tp,plotdat,banditcolz);
    
    % Significance testing
    if do.signif
        for b = 1:2
            p_crit          = .005;
            nit             = 1000;
            [p,praw]        = ClusterCorrection2(plotdat(:,:,b),nit,p_crit);
            stats.h(b,:)    = double(p <= p_crit);
            stats.h(stats.h == 0) = nan;
            stats.h(stats.h == 1) = 0;
            
            % Plot significance line
            plot(tp,stats.h(b,:)+sigpos(b),'s','MarkerFaceColor',banditcolz(b,:),'MarkerEdgeColor',banditcolz(b,:),'MarkerSize',6);
        end
    end
    
    % Legend
    hL = legend(steh,modlabels,'FontSize',labfntsz,'Location','NorthWest');
    hL.Position(1) = hL.Position(1)+.03;
    hL.Position(2) = hL.Position(2)+.04;
    legend boxoff
    
    % Additional settings and labels
    xlim([minX maxX]);
    ylim([minY maxY]);
    ax = gca;
    set(ax,'FontSize',16,'LineWidth',1.5);
    xlabel('Time from sample onset (ms)');
    ylabel('Dissimilarity correlation (Kendall \tau_A)');
    
    ax.XLabel.FontSize = labfntsz;
    ax.YLabel.FontSize = labfntsz;
end
