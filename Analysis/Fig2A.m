%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 2A: cross-validation RDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For Figure 1 - supplement 1C, change params.disttype to 'pearson'.

%% LOAD DATA

clc
clear

% Paths
savefolder          = 'Crossvalidation_RDM'; % folder to save newly created data
figfolder           = 'RSA_crossval'; % folder to save figures to
params.whichPhase   = 'test'; % use test phase data
params.disttype     = 'euclidean'; % Fig2A -> 'euclidean', Fig1-suppl1C -> 'pearson'

% Load stuff
Config_plot; % load plot variables

% Logicals
do.RDM              = true; % create RDM data, if it doesn't already exist
do.saveRDM          = true; % also save the RDM data, if it doesn't already exist

do.dimReduction     = false; % needs to be defined, but not used for this analysis
do.smooth           = true; % smooth the data?
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

%% Create cross-validation RDM

if do.RDM
    
    % Load Numbers data
    Numbers_load;
    num_data            = data;
    params.num_conds    = stim.samples; % Different conditions of RDM
    num_paths           = paths;
    
    % Load Bandit data
    Bandit_load;
    donk_data           = data;
    params.donk_conds   = subjranks; % Different conditions of RDM
    donk_paths          = paths;
    
    for s = 1:params.nsubj
        CreateCrossvalRDM(s,num_data,donk_data,do,params,num_paths,donk_paths);
    end
end

%% Plot cross-validation neural RDM

condnamez   = {'1','2','3','4','5','6';'b1','b2','b3','b4','b5','b6'};
eegrdm      = zeros(12,12,204,204,params.nsubj);

for s = 1:params.nsubj
    
    fprintf('\nConcatenate subject %d.\n',params.submat(s));
    
    % Load data
    inputfile = sprintf('Crossval_%03d_RDM_%s_%s',params.submat(s),params.whichPhase,params.disttype);
    load(fullfile(paths.data.save,inputfile));
    
    % Smoothing RDM
    if do.smooth
        wdwsz  	= 60/4;
        rdm.data = smoothRDM(rdm.data,wdwsz);
    end
    
    eegrdm(:,:,:,:,s) = rdm.data;
    
end

%% Plot RDMs

timewdw     = [350 600; 350 600];
numidx      = rdm.timepoints >= timewdw(1,1) & rdm.timepoints <= timewdw(1,2);
donkidx     = rdm.timepoints >= timewdw(2,1) & rdm.timepoints <= timewdw(2,2);
plotdat     = squeeze(mean(mean(mean(eegrdm(:,:,donkidx,numidx,:),3),4),5));

% Plot full 12 x 12 RDM
tmp         = sort(plotdat(:));
maplims     = [min(tmp(tmp > 0)) max(tmp)];
rdm.conds   = {'1','2','3','4','5','6','b1','b2','b3','b4','b5','b6'};

figB = figure;
colormap(hot);
imagesc(plotdat,maplims);
hc = colorbar;
axis square
ax = gca;
set(gca,'XTick',1:length(rdm.conds),'XTickLabel',rdm.conds,'FontSize',labfntsz);
set(gca,'YTick',1:length(rdm.conds),'YTickLabel',rdm.conds,'FontSize',labfntsz);
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ylabel(hc,'Correlation distance','FontSize',labfntsz);

% Plot cross-validation lower left quadrant
crossdat    = plotdat(7:12,1:6);
maplims     = [min(crossdat(:)) max(crossdat(:))];
rdm.conds   = {'1','2','3','4','5','6';'b1','b2','b3','b4','b5','b6'};

figB = figure;
colormap(hot);
imagesc(crossdat,maplims);
cax = colorbar('southoutside');
axis square
ax = gca;
set(ax,'XTick',1:length(rdm.conds),'XTickLabel',rdm.conds(1,:),'FontSize',labfntsz);
set(ax,'YTick',1:length(rdm.conds),'YTickLabel',rdm.conds(2,:),'FontSize',labfntsz);
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
xlabel(hc,'Euclidean distance','FontSize',axfntsz);
cax.FontSize = 24;







