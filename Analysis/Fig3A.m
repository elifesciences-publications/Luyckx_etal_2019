%%%%%%%%%%%%%%
%% Fig 3A: EDI
%%%%%%%%%%%%%%

%% LOAD DATA

clc
clear

% Paths
savefolder          = 'Crossvalidation_RDM'; % folder to save newly created data
figfolder           = 'RSA_crossval'; % folder to save figures to
params.whichPhase   = 'test'; % use test phase data
params.disttype     = 'euclidean';

% Load stuff
Config_plot; % load plot variables

% Logicals
do.RDM              = false; % create RDM data, if it doesn't already exist
do.saveRDM          = false; % also save the RDM data, if it doesn't already exist

do.dimReduction     = false; % needs to be defined, but not used for this analysis
do.smooth           = true; % smooth the data?

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

%% Cross-validation

clear models

% Obtain model RDMs
mods    = ModelRDM;
models  = fieldnames(mods);

% Subject correlations
for s = 1:params.nsubj
    
    fprintf('\nEDI lower left quadrant subject %d, %s phase.\n',params.submat(s),params.whichPhase);
    
    if s == 1
        edi = zeros(params.nsubj,204,204);
    end
    
    % Load data
    inputfile = sprintf('Crossval_%03d_RDM_%s_%s',params.submat(s),params.whichPhase,params.disttype);
    load(fullfile(paths.data.save,inputfile));
    
    % Smoothing RDM
    if do.smooth
        wdwsz    = 60/4; % size convolution kernel
        rdm.data = smoothRDM(rdm.data,wdwsz);
    end
    
    % Exemplar discriminability index (EDI)
    for t1 = 1:length(rdm.timepoints)
        for t2 = 1:length(rdm.timepoints)
            rdmvec = rdm.data(7:12,1:6,t1,t2);
            edi(s,t1,t2) = mean(rdmvec(~eye(6))) - mean(rdmvec(logical(eye(6)))); % off-diagonal - on-diagonal
        end
    end
    
end


%% Plot EDI

% Baseline correct
testdat     = edi;
tmp         = testdat;
tmp(:,rdm.timepoints >= 0, rdm.timepoints >= 0) = nan;
baseav      = nanmean(reshape(tmp,[params.nsubj,length(rdm.timepoints)^2]),2);

testdat     = testdat - baseav;

% Get t-values
[h,p,ci,stats]  = ttest(testdat);
tdat            = squeeze(stats.tstat);

mapval  = max(makeLong(abs(minmax(tdat))));
maplims = [linspace(-mapval,mapval,10)];
   
% Cluster test
nit         = 1000;
p_crit      = .005;
p_thresh    = .005;

fprintf('\nRunning cluster correction (n = %d), p_crit = %s, p_thresh = %s.\n',nit,num2str(p_crit),num2str(p_thresh));

[p,praw]    = ClusterCorrection2(testdat, nit, p_crit,p_thresh);
pmask       = double(squeeze(p <= p_thresh)); % 0 = not significant

% Plot figure
figC = figure;
colormap(hot);

[datah, ch] = contourf(rdm.timepoints,rdm.timepoints,tdat.*pmask,6); hold on
plot([0 0],[rdm.timepoints(1) rdm.timepoints(end)],'Color',[1 1 1]*.6,'LineWidth',lnwid);
plot([rdm.timepoints(1) rdm.timepoints(end)],[0 0],'Color',[1 1 1]*.6,'LineWidth',lnwid);
plot([rdm.timepoints(1) rdm.timepoints(end)],[rdm.timepoints(1) rdm.timepoints(end)],'Color',[1 1 1]*.6,'LineWidth',lnwid);

caxis([0 6]);
hcol = colorbar;
hcol.TickLabels{end} = '>6';
xlim([rdm.timepoints(1) rdm.timepoints(end)]);
ylim([rdm.timepoints(1) rdm.timepoints(end)]);

ax = gca;
axis square xy
set(ax,'FontSize',16,'LineWidth',1.5);
set(ax, 'Ticklength', [0 0]);

xlabel('Numerical: time (ms)','FontSize',labfntsz);
ylabel('Bandit: time (ms)','FontSize',labfntsz);
ylabel(hcol,'t-values','FontSize',labfntsz);
