%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 5C: SVD dimensionality reduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DATA

clc
clear

% Paths
savefolder          = 'Crossvalidation_RDM'; % folder to save newly created data
figfolder           = 'RSA_crossval'; % folder to save figures to
params.whichPhase   = 'test'; % use test phase data
params.disttype     = 'euclidean';
params.ndims        = 6; % number of dimensions


% Load stuff
Config_plot; % load plot variables

% Logicals
do.RDM              = true; % create RDM data, if it doesn't already exist
do.saveRDM          = true; % also save the RDM data, if it doesn't already exist

do.modelcorr        = true; %  
do.signif           = true; % run ANOVA on mean CV data
do.smooth           = true; % smooth the data?
do.plotting         = true; % plot?
do.save_data        = true; % save mean CV per dimensionality?

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
    
    % Combinations of subjects and dimensions
    allcombinations = combvec(1:params.nsubj,1:params.ndims)';
    
    for combo = 1:length(allcombinations)
        CreateCrossvalRDM_dim(allcombinations(combo,:),num_data,donk_data,do,params,num_paths,donk_paths);
    end
end

%% Cross-validation

% Obtain model RDMs
mods    = ModelRDM;
models  = fieldnames(mods);

if do.modelcorr
    for d = 1:params.ndims
        for s = 1:params.nsubj
            
            fprintf('\nCorrelating model - EEG RDM subject %d, %s phase, dim %d.\n',params.submat(s),params.whichPhase,d);
            
            if s == 1 && d == 1
                modcont = zeros(params.nsubj,204,204,params.ndims);
            end
            
            % Load data
            inputfile = sprintf('Crossval_%03d_RDM_%s_%s_dim%d',params.submat(s),params.whichPhase,params.disttype,d);
            load(fullfile(paths.data.saveEEG,inputfile));
            
            % Smoothing RDM
            if do.smooth
                wdwsz    = 60/4; % size convolution kernel
                rdm.data = smoothRDM(rdm.data,wdwsz);
            end
            
            % EEG - Model correlation
            actmod  = makeLong(mods.(models{1}));
            for t1 = 1:length(rdm.timepoints)
                for t2 = 1:length(rdm.timepoints)
                    rdmvec = makeLong(rdm.data(7:12,1:6,t1,t2));
                    modcont(s,t1,t2,d) = rankCorr_Kendall_taua(actmod,rdmvec);
                end
            end
        end
    end
    
    if do.save_data
        timepoints = rdm.timepoints;
        save(fullfile(paths.data.save,'Crossval_correlations_SVD'),'modcont','timepoints');
        fprintf('\nData saved.\n');
    end
    
else
    load(fullfile(paths.data.save,'Crossval_correlations_SVD'));
end

%% Plot all components of figure

if do.plotting
    
    %% Stats
    
    timewdw     = [350 600; 350 600];
    numidx      = timepoints >= timewdw(1,1) & timepoints <= timewdw(1,2);
    donkidx     = timepoints >= timewdw(2,1) & timepoints <= timewdw(2,2);
    
    for d = 1:params.ndims
        testdat(:,d) = squeeze(mean(mean(modcont(:,numidx,donkidx,d),2),3));
    end
    
    % Repeated-measures ANOVA
    if do.signif
        t = table([1:params.nsubj]',testdat(:,1),testdat(:,2),testdat(:,3),testdat(:,4),testdat(:,5),testdat(:,6));
        rm = fitrm(t,'Var2-Var7~Var1');
        ranova(rm)
    end
    
    % Pairwise comparison
    pairz = nchoosek(1:params.ndims,2);
    
    for tt = 1:length(pairz)
        [h(tt),p(tt),~,stats(tt)] = ttest(testdat(:,pairz(tt,1)),testdat(:,pairz(tt,2)));
    end
    
    %% P-value grid
    
    figP = figure;
    
    tdat = squareform(abs([stats.tstat]));
    
    colormap(cbrewer('seq','YlGn',100));
    maplims = [0 4];
    
    hm = imagesc(tdat,maplims);
    mask = 1-triu(ones(params.ndims,params.ndims));
    set(hm,'AlphaData',mask);
    ax = gca;
    box off
    
    plotvals = makeLong(squareform(p));
    textStrings = [repmat('p = ',36,1) num2str(plotvals, '%0.3f')];       % Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    textStrings(makeLong(logical(eye(params.ndims)))) = {''};
    textStrings(makeLong(logical(triu(ones(params.ndims,params.ndims))))) = {''};
    textStrings(strcmp(textStrings,'0.000')) = {'<0.001'};
    [x, y] = meshgrid(1:params.ndims);  % Create x and y coordinates for the strings
    hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
        'HorizontalAlignment', 'center');
    midValue = mean(get(ax, 'CLim'));  % Get the middle value of the color range
    textColors = repmat(makeLong(tdat) > midValue, 1, 3);  % Choose white or black for the text color of the strings so they can be easily seen over the background color
    set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors
    set(hStrings, 'FontSize', 18);
    set(hStrings, 'FontWeight', 'bold');

    set(ax, 'XTick',[],...
        'YTick', 2:6, ...                             % Change the axes tick marks
        'YTickLabel', {'2D','3D','4D','5D','6D'}, ...  %   and tick labels
        'TickLength', [0 0],...
        'FontSize',28);

    ylim([1.5 params.ndims+.5]);
    cax = colorbar;
    ylabel(cax,'t-values');
    set(cax,'Ticks',maplims(1):maplims(end),'TickLabels',{'0','1','2','3','>4'});
    
    %% Mean CV line plot (bottom)

    figL    = figure; hold on;
    ylims   = [.1 .18];
    meandat = mean(testdat);
    semdat  = std(testdat)./sqrt(params.nsubj);

    plot(1:6,meandat,'ks-','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','k');
    plot(1:6,meandat-semdat,'k-','LineWidth',2);
    plot(1:6,meandat+semdat,'k-','LineWidth',2);
    
    ax = gca;
    set(ax,'XTick',1:6,'XTickLabel',{'1D','2D','3D','4D','5D','6D'},'FontSize',28);
    set(ax,'YTick',linspace(ylims(1),ylims(end),5));
    ax.YAxis(1).FontSize = labfntsz;
    ylabel({'Mean CV'},'FontSize',28,'FontWeight','bold');
    
    xlim([.5 6.5]);
    ylim(ylims);
    
end
