%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 4CD: multiple regression w behavioural choice RDMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DATA

clc
clear

% Paths
savefolder          = 'Bandit_RDM'; % folder to save newly created data
figfolder           = 'Behaviour'; % folder to save figures to
params.whichPhase   = 'test'; % use test phase data
params.disttype     = 'euclidean';

% Load stuff
Config_plot; % load plot variables

% Logicals
do.modelcorr        = true; % correlate model RDMs with eeg?
do.smooth           = true; % smooth the data?
do.signif           = true;

%% Extra variables

% Load numbers data
Numbers_load;

% Numbers RDM (behaviour)
subresp = zeros(params.nsubj,params.nstim);
for s = 1:params.nsubj
    
    idx         = data.sub == params.submat(s) & data.r > 0;
    frame       = 1-unique(data.frame(idx));
    
    % Get response probabilities per number
    tmpSubResp  = makeLong(repmat(data.chosenCat(idx)-1,1,params.nsamp)');
    tmpCat      = makeLong(stim.category(idx,:)');
    tmpSamp     = makeLong(stim.samples(idx,:)');
    
    for c = 1:params.nstim
        subresp(s,c) = mean([tmpSubResp(tmpSamp == c & tmpCat == 2); -(tmpSubResp(tmpSamp == c & tmpCat == 1))+1]);
    end
    
end

subresp(params.frame == 0,:) = 1-subresp(params.frame == 0,:);
for s = 1:params.nsubj
    mods.num(:,:,s) = dist(subresp(s,:));
end

% Load bandit data
Bandit_load;

% Load RL fits
load(fullfile(paths.data.model,'Modelfit_full_test_RL'));

% Replace bandit index with perceived index
stim.comboSub = 0.*stim.combo;

for t = 1:params.ttrials
    for i = 1:params.nsamp
        stim.comboSub(t,i) = find(stim.combo(t,i) == mod.ranks(t,:));
    end
end

% Donkey behavioural RDM
highest = stim.comboSub(:,1) < stim.comboSub(:,2); % see whether subjective highest bandit is first (0) or second (1)

chosen_order     = 0*data.r - 99;
chosen_order(data.rmap == data.r) = 0; % chose first
chosen_order(data.rmap ~= data.r) = 1; % chose second
chosen_order(data.r < 0) = -99;

cor = highest == chosen_order;

mods.donkey = zeros(6,6,params.nsubj);
for s = 1:params.nsubj
    for i1 = 1:6
        for i2 = 1:6
            idx = data.sub == params.submat(s) & data.r > 0 & ((stim.comboSub(:,1) == i1 & stim.comboSub(:,2) == i2) | (stim.comboSub(:,1) == i2 & stim.comboSub(:,2) == i1));
            if sum(idx) == 0
                mods.donkey(i1,i2,s) = 0;
            else
                mods.donkey(i1,i2,s) = mean(cor(idx));
            end
        end
    end
end

%% Multiple regression with behavioural RDMs

if do.modelcorr
    
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
            elseif c == 2
                fprintf('\nBandit task: regression model - EEG RDM subject %d, %s phase.\n',params.submat(s),params.whichPhase);
                inputfile  = sprintf('Donkey_%03d_%s_EEG_RDM_%s',params.submat(s),params.whichPhase,params.disttype);
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
                vecs        = squeeze(allmod(:,:,s));
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
        modlabels   = {'Behavior numerical','Behavior bandit'};
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
        steh = myfillsteplot(rdm.timepoints(timeidx),plotdat,colz);
        
        % Significance testing
        if do.signif
            for m = 1:nmods
                p_crit          = .005;
                [p,praw]        = ClusterCorrection2(plotdat(:,:,m),5000,p_crit);
                stats.h(m,:)    = double(p <= p_crit);
                stats.h(stats.h == 0) = nan;
                stats.h(stats.h == 1) = 0;
                
                % Plot significance line
                plot(rdm.timepoints(timeidx),stats.h(m,:)+sigpos(m),'s','MarkerFaceColor',colz(m,:),'MarkerEdgeColor',colz(m,:),'MarkerSize',6);
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
end