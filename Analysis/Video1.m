%%%%%%%%%%%%%%%%
%% Video 1: MDS
%%%%%%%%%%%%%%%%

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
do.save_video       = false; % save video?

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

%% Make video

mds = nan(12,3,204);

% Get data on identical timepoints of both
for t = 1:length(rdm.timepoints)
    plotdat = squeeze(mean(mean(mean(mahalRSA(:,:,t,t,:),3),4),5));
    [mds(:,:,t), e] = cmdscale(flipfold(plotdat,1),3);
end

% Find minimum distance between adjacent configurations (for stability)
rot = [1 1; 1 -1; -1 1; -1 -1];

for t = 1:length(rdm.timepoints)-1
   
    for r = 1:4
        for i = 1:size(mds,1)
            distz(i) = pdist2(mds(i,1:2,t),rot(r,:).*mds(i,1:2,t+1));
        end
        mindist(r) = sum(distz);
    end
    
    mds(:,1:2,t+1) = rot(mindist == min(mindist),:).* mds(:,1:2,t+1);
end

% Plot frames for video
figure;
set(gcf, 'Position', [.1, .1, 600, 500]);
set(gcf,'color','w');

xlims = [-1 1];
ylims = [-1 1];

% create the video writer
if do.save_video
    v           = VideoWriter(fullfile(paths.figures.current,'MSD_1D.avi'));
    v.Quality   = 100;
    v.FrameRate = 16;

    open(v);
end

for t = 1:length(rdm.timepoints)
   
    hP(2) = plot(mds(7:12,1,t),mds(7:12,2,t),'o-','Color',colz(2,:),'MarkerFaceColor',colz(2,:),'LineWidth',lnwid,'MarkerSize',24);
    text(mds(7:12,1,t),mds(7:12,2,t),num2str([1:params.nstim]'),'Color','w','FontSize',22,'FontWeight','bold','HorizontalAlignment','center');   
    
    hold on;
    
    hP(1) = plot(mds(1:6,1,t),mds(1:6,2,t),'o-','Color',colz(1,:),'MarkerFaceColor',colz(1,:),'LineWidth',lnwid,'MarkerSize',24);
    text(mds(1:6,1,t),mds(1:6,2,t),num2str([1:params.nstim]'),'Color','w','FontSize',22,'FontWeight','bold','HorizontalAlignment','center');
        
    xlim(xlims);
    ylim(ylims);
    
    ax = gca;
    set(ax,'FontSize',axfntsz,'XTick',[],'YTick',[]);
    xlabel('Dim 1 (a.u.)','FontSize',labfntsz);
    ylabel('Dim 2 (a.u.)','FontSize',labfntsz);
    axis square
    
    hL = legend(hP,{'Numbers','Bandits'},'FontSize',labfntsz,'Location','NorthEast');
    hL.Position = [0.6300 0.7470 0.2050 0.1550];
    legend boxoff
    
    title(sprintf('%d ms',rdm.timepoints(t)),'FontSize',30);
        
    hold off
    drawnow
        
    % Save frame for video
    if do.save_video
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
end

% close the writer object
if do.save_video
    close(v);
end