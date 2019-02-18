%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 4AB: behavioural choice RDMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DATA

clc
clear

% Paths
savefolder          = 'Bandit_behaviour'; % folder to save newly created data
figfolder           = 'Behaviour'; % folder to save figures to
params.whichPhase   = 'test'; % use test phase data

% Load stuff
Config_plot; % load plot variables

% Logicals
do.signif   = false;

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

%% Behavioural correlations

if do.signif
    
    % Choice RDM after subtracting group mean
    num_resid       = mods.num - repmat(mean(mods.num,3),1,1,params.nsubj);
    donkey_resid    = mods.donkey - repmat(mean(mods.donkey,3),1,1,params.nsubj);
    
    for s = 1:params.nsubj
        corr_full(s) = rankCorr_Kendall_taua(squareform(mods.num(:,:,s)),squareform(mods.donkey(:,:,s)));
        corr_resid(s) = rankCorr_Kendall_taua(squareform(num_resid(:,:,s)),squareform(donkey_resid(:,:,s)));
    end
    
    [p,h,stats] = signrank(corr_full)
    [p,h,stats] = signrank(corr_resid)
end

%% Plot behavioural RDMs

% Numerical RDM
figN = figure;
colormap(hot);
plotdat = mean(mods.num,3);
mapvals = [min(plotdat(plotdat ~= 0)) max(plotdat(:))];
imagesc(plotdat,mapvals);

ax  = gca;
axis square
set(gca,'XTick',1:6,'XTickLabel',num2str([1:6]'));
set(gca,'YTick',1:6,'YTickLabel',num2str([1:6]'));
ax.XAxis.FontSize   = 24;
ax.YAxis.FontSize   = 24;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';

cax = colorbar;
set(cax,'FontSize',20);
ylabel(cax,'\Delta P(choice)','FontSize',24);

% Bandit RDM
figD = figure;
colormap(hot);
plotdat = mean(mods.donkey,3);
mapvals = [min(plotdat(plotdat ~= 0)) max(plotdat(:))];
imagesc(plotdat,[.73 .95]);

ax  = gca;
axis square
set(gca,'XTick',1:6,'XTickLabel',{'b1','b2','b3','b4','b5','b6'});
set(gca,'YTick',1:6,'YTickLabel',{'b1','b2','b3','b4','b5','b6'});
ax.XAxis.FontSize   = 24;
ax.YAxis.FontSize   = 24;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';

cax = colorbar;
set(cax,'FontSize',20);
ylabel(cax,'P(highest)','FontSize',24);