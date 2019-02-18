%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 1 - suppl 2A-B: RSA individual tasks - split per framing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DATA

clc
clear

% Paths
savefolder          = 'Bandit_RDM';
figfolder           = 'RSA_indiv'; % folder to save figures to
params.whichPhase   = 'test'; % use test phase data
params.disttype     = 'euclidean';

% Load stuff
Config_plot; % load plot variables

% Logicals
do.RDM              = false; % create RDM data, if it doesn't already exist
do.saveRDM          = false; % also save the RDM data, if it doesn't already exist

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

%% Define frames

whichFrame  = params.frame; % 1 = high frame, 0 = low frame
framenames  = {'high','low'};
frameid     = [1 0];
nframes     = length(unique(whichFrame));

%% Plot data

if do.plotting
    
    %% RSA lines
    
    % Settings
    nmods       = size(plotdat,3);
    modlabels   = {'Numerical','Bandit'};
    sigtimes    = zeros(2,ntimepoints,nframes);
    
    ylims       = [-.1 .4];
    sigpos      = [-.08, -.07]; % Y position of significance
    
    minX        = min(mutual_time);
    maxX        = max(mutual_time);
    minY        = ylims(1);
    maxY        = ylims(2);
       
    for f = 1:nframes % loop per task framing or collapsed
        
        % Subset of participants per task frame
        fidx   	= whichFrame == frameid(f);

        % Figure for line plot
        figR(f) = figure; hold on;
        
        % Plot lines
        plot([minX maxX],[0 0],'k--','LineWidth',1.5);
        plot([0 0],[minY*1.1 maxY*1.1],'k-','LineWidth',1.5);
        steh = myfillsteplot(mutual_time,plotdat(fidx,:,:),colz);
        
        % Significance testing
        if do.signif
            for m = 1:nmods
                p_crit          = .005;
                nit             = 1000;
                [p,praw]        = ClusterCorrection2(plotdat(fidx,:,m),nit,p_crit);
                sigtimes(m,:,f) = double(p <= p_crit);
                stats.h(m,:)    = sigtimes(m,:,f);
                stats.h(m,stats.h(m,:) == 0) = nan;
                stats.h(m,stats.h(m,:) == 1) = 0;
                
                % Plot significance line
                plot(mutual_time,stats.h(m,:)+sigpos(m),'s','MarkerFaceColor',colz(m,:),'MarkerEdgeColor',colz(m,:),'MarkerSize',6);
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
        
        title(framenames{f});
    end
    

end
