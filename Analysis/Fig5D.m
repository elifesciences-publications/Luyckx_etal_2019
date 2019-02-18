%%%%%%%%%%%%%%
%% Fig 5D: MDS
%%%%%%%%%%%%%%

%% LOAD DATA

clc
clear

% Paths
savefolder          = 'Crossvalidation_RDM'; % folder to save newly created data
figfolder           = 'MDS'; % folder to save figures to
params.whichPhase   = 'test'; % use test phase data bandit task
params.disttype     = 'euclidean';

% Load stuff
Config_plot; % load plot variables

% Logicals
do.RDM              = false; % create RDM data, if it doesn't already exist
do.saveRDM          = false; % also save the RDM data, if it doesn't already exist

do.smooth           = true; % smooth data?

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

%% Load RDM data

mahalRSA = zeros(12,12,204,204,params.nsubj);

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
    
    mahalRSA(:,:,:,:,s) = rdm.data;
    
end

%% Plot

timewdw     = [350 600; 350 600];
numidx      = rdm.timepoints >= timewdw(1,1) & rdm.timepoints <= timewdw(1,2);
donkidx     = rdm.timepoints >= timewdw(2,1) & rdm.timepoints <= timewdw(2,2);
plotdat     = squeeze(mean(mean(mean(mahalRSA(:,:,numidx,donkidx,:),3),4),5));
plotdat(logical(eye(12))) = 0; % zero the diagonal to allow for mds

xlims = [-.6 1];
ylims = [-.5 .8];

figD = figure; hold on;
[mds e] = cmdscale(plotdat);

hP(2) = plot(mds(7:12,1),mds(7:12,2),'o-','Color',colz(2,:),'MarkerFaceColor',colz(2,:),'LineWidth',lnwid,'MarkerSize',20);
text(mds(7:12,1),mds(7:12,2),num2str([1:params.nstim]'),'Color','w','FontSize',16,'FontWeight','normal','HorizontalAlignment','center');
hP(1) = plot(mds(1:6,1),mds(1:6,2),'o-','Color',colz(1,:),'MarkerFaceColor',colz(1,:),'LineWidth',lnwid,'MarkerSize',20);
text(mds(1:6,1),mds(1:6,2),num2str([1:params.nstim]'),'Color','w','FontSize',16,'FontWeight','normal','HorizontalAlignment','center');

xlim(xlims);
ylim(ylims);

ax = gca;
set(ax,'FontSize',axfntsz,'XTick',[],'YTick',[]);
xlabel('Dim 1 (a.u.)','FontSize',labfntsz);
ylabel('Dim 2 (a.u.)','FontSize',labfntsz);
axis square

hL = legend(hP,{'Numbers','Bandits'},'FontSize',labfntsz,'Location','NorthEast');
legend boxoff
box
