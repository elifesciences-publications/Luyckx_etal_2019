%%%%%%%%%%%%%%%%%
%% Fig 5 A-B: CPP
%%%%%%%%%%%%%%%%%

%% LOAD DATA

clc
clear

% Paths
savefolder          = ''; % folder to save newly created data
figfolder           = 'CPP'; % folder to save figures to
params.whichPhase   = 'test'; % use test phase data bandit task

% Load stuff
Config_plot; % load plot variables

% Logicals
do.plotting         = false; % plot?
do.save_bar         = false; % save peak data?

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

% Get mutual time
mutual_time = [-64 748];

%% Concatenate data

for c = 1:2
    
    for s = 1:params.nsubj
        
        fprintf('\n%s: processing ERP subject %d.\n',params.tasknames{c},params.submat(s));
        
        if c == 1
            
            % Load paths/data
            if s == 1
                Numbers_load;
            end
            
            % File name
            inputfile   = sprintf('Numbers_%03d_',params.submat(s));
            
            % Different conditions
            data.conds  = stim.samples;
            
        elseif c == 2
            
            % Load paths/data
            if s == 1
                Bandit_load;
            end
            
            % File name
            inputfile   = sprintf('Donkey_%03d_%s_',params.submat(s),params.whichPhase);
            
            % Different conditions
            data.conds  = subjranks;
        end
        
        %% Initialise
        
        % Index for exluding trials
        goodtrials  = [];
        bindx       = []; % behavioural index
        eindx       = []; % eeg index
        
        %% Load data
        
        load(fullfile(paths.data.EEG.raw,sprintf('%sEEG_samples.mat',inputfile)));
        load(fullfile(paths.data.EEG.raw,sprintf('%srejectedTrials.mat',inputfile)));
        
        %% Good trials
        
        goodtrials  = 1-rejectedTrialz;
        
        % Variables
        ntimepoints = length(eeg.timepoints); % number of time points in epoch
        rtrials     = length(find(goodtrials == 1)); % number of preserved trials
        
        if s == 1
            sampdat = zeros(params.nsubj,eeg.nbchan,length(eeg.timepoints),length(unique(data.conds)));
        end
        
        %% Indices
        
        idx     = data.sub == params.submat(s);
        cortmp  = data.RT(idx) > 0;
        bindx   = cortmp & goodtrials'; % behavioural index
        
        eindx   = cortmp(logical(goodtrials')) == 1; % eeg index
        eindx   = makeLong(repmat(eindx,1,params.nsamp)');
        
        %% Averaged ERP
        
        tmpsamp = data.conds(idx,:);
        samps   = makeLong(tmpsamp(bindx,:)');
        tmpY    = zscore(eeg.data(:,:,eindx),[],3);
        
        for n = 1:params.nstim
            sampdat(s,:,:,n) = squeeze(mean(tmpY(:,:,samps == n),3));
        end
        
    end
    
    % Save normalised data with right timing
    timeidx             = eeg.timepoints >= mutual_time(1) & eeg.timepoints <= mutual_time(end);
    normdat(:,:,:,:,c)  = sampdat(:,:,timeidx,:);
    
end

%% Get CPP data

cpp     = {'CP1','P1','POZ','PZ','CPZ','CP2','P2'};
chanidx = label2index(cpp,{chanlocs.labels}); % Get indices of electrodes to plot

cppdat  = squeeze(mean(normdat(:,chanidx,:,:,:),2));

%% Peak amplitude per sample

% Identify significant time window (leave-one-out)
p_crit      = .01; % critical p-value for FDR correction
myclust     = 50; % min size of cluster
timewdw     = zeros(params.nsubj,2,2);
bardat      = zeros(params.nsubj,6,2);
timepoints  = mutual_time(1):4:mutual_time(end);
xnamez      = {'Digit','Bandit (worst to best)'};

for c = 1:2 % loop over tasks
    for s = 1:params.nsubj
        
        submat  = setdiff(1:params.nsubj,s); % leave one out
        tmpcpp  = cppdat(:,:,:,c);
        subdat 	= tmpcpp(submat,:,:); % get remaining subjects data
        
        % Omnibus test per timepoint
        for t = 1:length(timepoints)
            testdat             = squeeze(subdat(:,t,:));
            [p(t),tbl,stats]    = kruskalwallis(testdat,[],'off');
        end
        
        % Cluster correction (find largest modulation)
        [p_fdr, p_masked]   = fdr(p,p_crit);
        supfin              = find(p_masked==1);
        p_clustmasked       = zeros(size(p_masked));
        
        for k = 2:length(supfin) %  cluster find
            if supfin(k) + myclust < length(timepoints)
                if  p_masked(supfin(k)) == 1 && sum(p_masked(supfin(k):supfin(k)+myclust-1)) == myclust
                    p_clustmasked(supfin(k):supfin(k)+myclust-1) = 1;
                end
            end
        end
        
        timez        	= find(p_clustmasked == 1);
        timewdw(s,1,c)	= timez(1);
        timewdw(s,2,c)  = timez(end);
        
        bardat(s,:,c)   = squeeze(mean(tmpcpp(s,timez,:)));
    end
    
    % Save peak data
    if do.save_bar
        
        % Load paths/data
        if c == 1
            Numbers_load;
        elseif c == 2
            Bandit_load;
        end
        
        % Select data
        peakz   = squeeze(bardat(:,:,c));
        
        save(fullfile(paths.data.EEG.raw,sprintf('Peak_CPP_%s_%s.mat',params.tasknames{c},'collapsed')),'peakz');
        fprintf('\nSaved CCP peak data %s, %s.\n',params.tasknames{c},'collapsed');
    end
end

%% Topoplot sample X at significant time window

% Open eeglab to reach eggheadplot plugin
eeglab;
close

whichSamp = 6; % which sample to plot (1-6)

for c = 1:2 % iterate over tasks
    whichTime   = round(median(timewdw(:,1,c))):round(median(timewdw(:,2,c))); % get activity from median time window
    plotdat     = squeeze(mean(mean(normdat(:,:,whichTime,whichSamp,c),1),3));
    maplims     = [-.1 .1]; % set scale
    chanidx     = label2index(cpp,{chanlocs.labels}); % Get indices of electrodes to plot
    
    % Egghead plot
    eggheadplot('Channels', {eeg.chanlocs.labels}, 'Amplitude', plotdat, 'Method', 'natural', 'Scale', maplims, 'Contours', 0, 'FillColor', [1 1 1], 'MapStyle', 'jet', 'Style', 'Full', 'ElectrodeSize', 10, 'ShowBrain','No','Smooth',50);
    axis equal
end

%% Plot ERP lines

leglabels   = {'6','5','4','3','2','1';'b6','b5','b4','b3','b2','b1'};

for c = 1:2
    ylims       = [-.12 .2];
    masstest    = 'cluster';
    signif      = 0;
    
    figR = figure;
    [steh] = fEEG_steplot(cppdat(:,:,:,c),timepoints,ylims,signif,masstest,erpcol,'bottom');
    
    minX = timepoints(ceil(mean(timewdw(:,1,c))));
    maxX = timepoints(ceil(mean(timewdw(:,2,c))));
    sigfill = fill([minX maxX maxX minX],[min(ylims),min(ylims) max(ylims) max(ylims)],[.6 .6 .6],'Edgecolor','none','FaceAlpha',.25);
    uistack(sigfill,'bottom');
    
    % Custom legend
    erpcol2     = flipud(erpcol);
    coord       = [700, .18];
    stepsz      = .02;
    for l = 1:length(leglabels)
        text(coord(1),coord(2)-(stepsz*(l-1)),leglabels{c,l},'Color',erpcol2(l,:),'FontSize',18,'FontWeight','bold');
    end
    
    set(gca,'FontSize',axfntsz);
    xlabel('Time from sample onset (ms)','FontSize',labfntsz);
    ylabel('Norm. amplitude','FontSize',labfntsz);
    
end
