%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 1B: RSA magnitude individual tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For Figure 1 - supplement 1A, change params.disttype to 'pearson'.

%% LOAD DATA

clc
clear

% Paths
savefolder          = '';
figfolder           = 'RSA_indiv'; % folder to save figures to
params.whichPhase   = 'test'; % use test phase data
params.disttype     = 'euclidean'; % Fig1B -> 'euclidean', Fig1-suppl1A -> 'pearson'

% Load stuff
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
Bandit_load;
load(fullfile(paths.data.model,'Modelfit_full_test_RL'));

% Replace bandit index with perceived index
subjranks = 0.*stim.combo;

for t = 1:params.ttrials
    for i = 1:params.nsamp
        subjranks(t,i) = find(stim.combo(t,i) == mod.ranks(t,:));
    end
end

% Labels
params.tasknames    = {'numbers','bandit'};
condnamez           = {'1','2','3','4','5','6';'b1','b2','b3','b4','b5','b6'};

%% Create RDMs

if do.RDM
    
    for c = 1:2 % loop through tasks
        for s = 1:params.nsubj
            
            if c == 1
                
                fprintf('\nCreating RDMs numerical task.\n');
                
                % Load paths/data
                if s == 1
                    Numbers_load;
                end
                
                % File names
                inputfile   = sprintf('Numbers_%03d_',params.submat(s));
                outputfile  = sprintf('Numbers_%03d_EEG_RDM_%s',params.submat(s),params.disttype);
                
                % Different conditions of RDM
                data.conds  = stim.samples;
                
            elseif c == 2
                
                fprintf('\nCreating RDMs bandit task.\n');
                
                % Load paths/data
                if s == 1
                    Bandit_load;
                end
                
                % File names
                inputfile  	= sprintf('Donkey_%03d_%s_',params.submat(s),params.whichPhase);
                outputfile  = sprintf('Donkey_%03d_%s_EEG_RDM_%s',params.submat(s),params.whichPhase,params.disttype);
                
                data.conds  = subjranks;
            end
            
            % Create RDMs
            CreateRDM(s,inputfile,outputfile,do,data,params,paths.data.EEG);
        end
    end
end

%% Create model RDM

mods 	= ModelRDM;
models  = fieldnames(mods);

%% RSA

if do.modelcorr
    
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
            
            % Initialise
            if c == 1 && s == 1
                rsa_results = zeros(params.nsubj,length(rdm.timepoints),2);
            end
            
            % EEG - Model correlation
            actmod = squareform(mods.(models{1}));
            for t = 1:length(rdm.timepoints)
                rsa_results(s,t,c) = rankCorr_Kendall_taua(actmod,squareform(rdm.data(:,:,t)));
            end
            
        end
        
        % Store time points
        if c == 1
            timepoints1 = rdm.timepoints;
        elseif c == 2
            timepoints2 = rdm.timepoints;
        end
    end
end

%% Get mutual time

mutual_time     = intersect(timepoints1,timepoints2);
ntimepoints     = length(mutual_time);
eeg1wdw         = timepoints1 >= mutual_time(1) & timepoints1 <= mutual_time(end); % numbers idx
eeg2wdw         = timepoints2 >= mutual_time(1) & timepoints2 <= mutual_time(end); % bandit idx

% Get data for line plots
plotdat         = zeros(params.nsubj,ntimepoints,2);
plotdat(:,:,1)  = rsa_results(:,eeg1wdw,1);
plotdat(:,:,2)  = rsa_results(:,eeg2wdw,2);

%% Plot data

if do.plotting
    
    %% RSA lines
    
    % Settings
    nmods       = size(plotdat,3);
    modlabels   = {'Numerical','Bandit'};
    sigtimes    = zeros(2,ntimepoints);
    
    ylims       = [-.07 .35];
    sigpos      = [-.04, -.05]; % Y position of significance
    
    minX        = min(mutual_time);
    maxX        = max(mutual_time);
    minY        = ylims(1);
    maxY        = ylims(2);
    
    % Figure for line plot
    figR = figure; hold on;
    
    % Plot lines
    plot([minX maxX],[0 0],'k--','LineWidth',1.5);
    plot([0 0],[minY*1.1 maxY*1.1],'k-','LineWidth',1.5);
    steh = myfillsteplot(mutual_time,plotdat,colz);
    
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
            plot(mutual_time,stats.h(b,:)+sigpos(b),'s','MarkerFaceColor',colz(b,:),'MarkerEdgeColor',colz(b,:),'MarkerSize',6);
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
    ylabel('Dissimilarity correlation (Kendall \tau-a)');
    
    ax.XLabel.FontSize = labfntsz;
    ax.YLabel.FontSize = labfntsz;
    
    %% Magnitude model
    
    figR = figure;
    colormap(hot);
    imagesc(mods.numd);
    set(gca,'XTick',[],'YTick',[]);
    axis square
    
end
