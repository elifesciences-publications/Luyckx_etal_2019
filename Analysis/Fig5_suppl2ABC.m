%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 5 - suppl2 A-C: CPP
%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Labels
condnamez = {'1','2','3','4','5','6';'b1','b2','b3','b4','b5','b6'};

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

%% Plot bar data

colmap  = erpcol;
ylims   = [-.08 .12];
inc     = .04;

for c = 1:2
    
    M       = mean(bardat(:,:,c));
    sembar  = squeeze(std(bardat(:,:,c))./sqrt(params.nsubj));
    
    figure; hold on;
    plot([0 params.nstim+1],[0 0],'k','LineWidth',1.5);
    
    for i = 1:params.nstim
        bar(i,M(i),.5,'FaceColor',colmap(i,:),'EdgeColor','k','LineWidth',1.5);
        errorbar(i,M(i),sembar(i),'k','MarkerSize',0.01,'LineWidth',1.5);
    end
    
    ax = gca;
    set(ax,'FontSize',axfntsz);
    ax.XAxis.FontSize   = labfntsz;
    ax.XAxis.FontWeight = 'bold';
    set(ax,'XTick',[1:params.nstim],'XTickLabel',condnamez(c,:));
    set(ax,'YTick',min(ylims):inc:max(ylims));
    xlabel(xnamez{c},'FontSize',labfntsz,'FontWeight','normal');
    ylabel('Norm. amplitude','FontSize',labfntsz);
    xlim([0 params.nstim+1]);
    ylim(ylims);

end

%% Plot CPP RDMs within task

for c = 1:2    
    
    for s = 1:params.nsubj
        rdm(:,:,s) = dist(bardat(s,:,c));
    end
    
    figure;
    plotdat     = squeeze(mean(rdm,3));
    tmp         = sort(plotdat(:));
    maplims     = [min(tmp(tmp > 0)) max(tmp)];
    colormap(hot);
    imagesc(plotdat,maplims);
    axis square;
    ax = gca;
    set(ax,'FontSize',labfntsz,'FontWeight','bold');
    set(ax,'XTick',1:6,'XTickLabel',condnamez(c,:));
    set(ax,'YTick',1:6,'YTickLabel',condnamez(c,:));
end

%% Plot CPP RDM between tasks
    
rdm_all = zeros(12,12,params.nsubj);

for s = 1:params.nsubj
    rdm_all(:,:,s) = dist(squeeze([bardat(s,:,1),bardat(s,:,2)]));
end

figure;
plotdat     = squeeze(mean(rdm_all,3));
tmp         = sort(plotdat(:));
maplims     = [min(tmp(tmp > 0)) max(tmp)];
colormap(hot);
imagesc(plotdat,maplims);

axis square;
ax = gca;
set(gca,'XTick',1:length(makeLong(condnamez)),'XTickLabel',makeLong(condnamez')');
set(gca,'YTick',1:length(makeLong(condnamez)),'YTickLabel',makeLong(condnamez')');
set(ax,'FontSize',lgndfntsz,'FontWeight','bold');

