%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 5 - suppl 2 DE: multiple regression control CPP RDMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DATA

clc
clear

% Paths
savefolder          = ''; % folder to save newly created data
figfolder           = 'CPP'; % folder to save figures to
params.whichPhase   = 'test'; % use test phase data
params.disttype     = 'euclidean';

% Load stuff
Config_plot; % load plot variables

% Logicals
do.RDM              = false; % create RDM data, if it doesn't already exist
do.saveRDM          = false; % also save the RDM data, if it doesn't already exist

do.smooth           = true; % smooth the data?
do.signif           = true; % run non-parametric cluster test

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
params.tasknames    = {'Numerical','Bandit'};
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

%% Multiple regression with CPP RDMs

% Loop over tasks
for c = 1:2
    
    clear modcont models
    
    % Load correct paths
    if c == 1
        Numbers_load;
    elseif c == 2
        Bandit_load;
    end
    
    % Subject correlations
    for s = 1:params.nsubj

        % Load subject data
        if c == 1
            fprintf('\nNumerical task: regression model - EEG RDM subject %d.\n',params.submat(s));
            inputfile  = sprintf('Numbers_%03d_EEG_RDM_%s',params.submat(s),params.disttype);
            
            % Obtain model RDMs
            mods    = ModelRDM_CPP(s,params);
            mods    = rmfield(mods,{'crossCPP','donkeyCPP'});
            models  = fieldnames(mods);
            
        elseif c == 2
            fprintf('\nBandit task: regression model - EEG RDM subject %d, %s phase.\n',params.submat(s),params.whichPhase);
            inputfile  = sprintf('Donkey_%03d_%s_EEG_RDM_%s',params.submat(s),params.whichPhase,params.disttype);
            
            % Obtain model RDMs
            mods    = ModelRDM_CPP(s,params);
            mods    = rmfield(mods,{'numCPP','crossCPP'});
            models  = fieldnames(mods);            
        end
        
        load(fullfile(paths.data.EEG.RDM,inputfile));

        % Smooth data
        if do.smooth
            wdwsz    = 60/4; % size convolution kernel
            rdm.data = smoothRDM(rdm.data,wdwsz);
        end
        
        % Regression (Pearson)
        models  = fieldnames(mods);
        for m = 1:length(models)
            allmod      = mods.(models{m});
            vecs        = squeeze(allmod);
            actmod(:,m) = zscore(squareform(vecs));
        end
        
        for t = 1:length(rdm.timepoints)
            betas = regress(zscore(squareform(rdm.data(:,:,t)))',[actmod(:,1)*0+1 actmod]);
            modcont(s,t,:) = betas(2:end);
        end
        
    end
    
    %% Get mutual time
    
    mutual_time = [-64 748];
    timeidx     = rdm.timepoints >= mutual_time(1) & rdm.timepoints <= mutual_time(end);
    
    %% Plot correlations
    
    
    % Settings
    ylims       = [-.2 .6];
    plotdat     = modcont(:,timeidx,:);
    semdat      = squeeze(nansem(plotdat));
    nmods       = size(plotdat,3);
    modlabels   = {params.tasknames{c},'CPP'};
    colmat      = [1,4; 2,4];
    sigpos      = [-.15, -.17]; % Y position of significance
    labfntsz    = 20;
    
    minX        = min(rdm.timepoints(timeidx));
    maxX        = max(rdm.timepoints(timeidx));
    minY        = ylims(1);
    maxY        = ylims(2);
    
    % Figure
    figR = figure; hold on;
    
    % Plots
    plot([minX maxX],[0 0],'k--','LineWidth',1.5);
    plot([0 0],[minY*1.1 maxY*1.1],'k-','LineWidth',1.5);
    steh = myfillsteplot(rdm.timepoints(timeidx),plotdat,colz(colmat(c,:),:));
    
    % Significance testing
    if do.signif
        for m = 1:nmods
            p_crit          = .005;
            [p,praw]        = ClusterCorrection2(plotdat(:,:,m),5000,p_crit);
            stats.h(m,:)    = double(p <= p_crit);
            stats.h(stats.h == 0) = nan;
            stats.h(stats.h == 1) = 0;
            
            % Plot significance line
            plot(rdm.timepoints(timeidx),stats.h(m,:)+sigpos(m),'s','MarkerFaceColor',colz(colmat(c,m),:),'MarkerEdgeColor',colz(colmat(c,m),:),'MarkerSize',6);
        end
    end
    
    % Legend
    hL = legend(steh,modlabels,'FontSize',labfntsz,'Location','Best');
    legend boxoff
    
    % Additional settings and labels
    xlim([minX maxX]);
    ylim([minY maxY]);
    ax = gca;
    set(ax,'FontSize',16,'LineWidth',1.5);
    xlabel('Time from sample onset (ms)');
    ylabel('Mean \beta-weight');
    
    ax.XLabel.FontSize = labfntsz;
    ax.YLabel.FontSize = labfntsz;
    
end