%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 2 - suppl 2: RSA cross-validation permutations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DATA

clc
clear

% Paths
savefolder          = 'Crossvalidation_RDM'; % folder to save newly created data
figfolder           = 'RSA_crossval'; % folder to save figures to
params.whichPhase   = 'test'; % use test phase data
params.disttype     = 'euclidean'; % Fig2B -> 'euclidean', Fig1-suppl1D -> 'pearson'

% Load stuff
Config_plot; % load plot variables

% Logicals
do.parpooling       = false; % when using supercomputer with multiple workers
do.RDM              = true; % create RDM data, if it doesn't already exist
do.saveRDM          = true; % also save the RDM data, if it doesn't already exist
do.zval             = true; % calculate z-values
do.save_zval        = true; % save z-values

do.smooth           = true; % smooth the data?
do.plotting         = true; % plot?

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

%% Create permutations RDMs (Numbers subject X -> Bandit subject Y)

if do.RDM
    try
        tic
        
        % Parpooling
        if do.parpooling
            numWorkers = 32; % memory too small to do 46 at same time
            c = parcluster();
            c.NumWorkers = numWorkers;
            parpool(c,numWorkers);
        else
            numWorkers = 0;
        end
        
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
        
        % Combinations of subjects
        allcombinations = combvec(1:params.nsubj,1:params.nsubj-1)';
        
        % Create RDM
        parfor (combo = 1:length(allcombinations),numWorkers)
            CreateCrossvalRDM_perm(allcombinations(combo,:),num_data,donk_data,do,params,num_paths,donk_paths);
        end
        toc
        
        % End parpool session
        if do.parpooling
            delete(gcp);
        end
        
    catch ME
        if do.parpooling
            delete(gcp);
        end
        rethrow(ME)
        toc
        return
    end
end

%% Get Z-values (distribution of correlations <> within-subject crossval)

if do.zval
    
    try
        tic
        
        % Parpooling
        if do.parpooling
            numWorkers = params.nsubj;
            c = parcluster();
            c.NumWorkers = numWorkers;
            parpool(c,numWorkers);
        else
            numWorkers = 0;
        end
        
        % Calculate z-values
        parfor (s = 1:params.nsubj,numWorkers)
            ZValues_perm(s,do,params,paths);
        end
        toc
        
        % End parpool session
        if do.parpooling
            delete(gcp);
        end
        
    catch ME
        if do.parpooling
            delete(gcp);
        end
        rethrow(ME)
        toc
        return
    end
end

%% Concatenate data

timepoints  = -64:4:748;
allZ        = zeros(params.nsubj,length(timepoints),length(timepoints));

for s = 1:params.nsubj
    
    fprintf('\nLoading Z-values sub %d\n',params.submat(s));
    
    load(fullfile(paths.data.save,sprintf('Zval_%03d_perm',params.submat(s))));
    
    allZ(s,:,:) = Z;
end

%% Plot correlations

if do.plotting
    
    testdat         = allZ;
    [h,p,ci,stats]  = ttest(testdat);
    tdat            = squeeze(stats.tstat);
    
    mapval  = max(makeLong(abs(minmax(tdat))));
    maplims = linspace(-mapval,mapval,10);
    maxT    = 4;
    
    % Cluster test
    nit         = 1000;
    p_crit      = .05;
    p_thresh    = .05;
    
    fprintf('\nRunning cluster correction (n = %d), p_crit = %s, p_thresh = %s.\n',nit,num2str(p_crit),num2str(p_thresh));
    
    [p,praw]    = ClusterCorrection2(testdat, nit, p_crit,p_thresh);
    pmask       = double(squeeze(p <= p_thresh)); % 0 is not significant
    
    % Plot figure
    figC = figure;
    colormap(hot);
    
    [datah, ch] = contourf(timepoints,timepoints,tdat.*pmask,4); hold on
    plot([0 0],[timepoints(1) timepoints(end)],'Color',[1 1 1]*.6,'LineWidth',lnwid);
    plot([timepoints(1) timepoints(end)],[0 0],'Color',[1 1 1]*.6,'LineWidth',lnwid);
    plot([timepoints(1) timepoints(end)],[timepoints(1) timepoints(end)],'Color',[1 1 1]*.6,'LineWidth',lnwid);
    
    caxis([0 maxT]);
    hcol = colorbar;
    hcol.TickLabels{end} = sprintf('>%d',maxT);
    xlim([timepoints(1) timepoints(end)]);
    ylim([timepoints(1) timepoints(end)]);
    
    ax = gca;
    axis square xy
    set(ax,'FontSize',16,'LineWidth',1.5);
    set(ax, 'Ticklength', [0 0]);
    
    xlabel('Numerical: time (ms)','FontSize',labfntsz);
    ylabel('Bandit: time (ms)','FontSize',labfntsz);
    ylabel(hcol,'t-values','FontSize',labfntsz);
    
end

