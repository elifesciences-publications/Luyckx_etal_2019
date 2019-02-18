%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 4EF: anti-compression of behaviour and EEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DATA

clear

% Paths
savefolder          = ''; % folder to save newly created data (unused)
figfolder           = 'Model';
params.whichPhase   = 'test'; % use test phase data
params.disttype     = 'euclidean'; % type of distance measure

% Load stuff
Config_plot; % load plot variables

% Logicals
do.smooth           = true; % smooth the data?
do.plotting         = true; % plot?
do.signif           = true; % test significance?
do.save_plot        = false; % save plot?

%% Extra variables

neurom 	= @(k,f) (abs(f).^k);
kspace  = .35:.01:3;
f       = linspace(0,1,6);
timewdw = [350 600];

expnamez = {'Numerical','Bandit'};

% Load behavioural data
Numbers_load;
load(fullfile(paths.data.model,sprintf('Numbers_modelfit')));
bestk.behav  = bestfit(:,2);

%% Find best k EEG

corrmap = zeros(length(kspace),params.nsubj,2);

for c = 1:2 % loop over tasks
    
    % Subject correlations
    for s = 1:params.nsubj
        
        fprintf('\nFind best k for EEG RDM subject %d, %s.\n',params.submat(s),expnamez{c});
        
        if c == 1
            % Load paths
            if s == 1
                Numbers_load;
            end
            inputfile   = sprintf('Numbers_%03d_EEG_RDM_%s',params.submat(s),params.disttype);
        elseif c == 2
            % Load paths
            if s == 1
                Bandit_load;
            end
            inputfile   = sprintf('Donkey_%03d_%s_EEG_RDM_%s',params.submat(s),params.whichPhase,params.disttype);
        end
        
        % Load data
        load(fullfile(paths.data.EEG.RDM,inputfile));
        
        if do.smooth
            wdwsz = 60/4;
            rdm.data = smoothRDM(rdm.data,wdwsz);
        end
        
        timez           = rdm.timepoints >= timewdw(1) & rdm.timepoints <= timewdw(2);
        eegrdm(:,:,s)	= mean(rdm.data(:,:,timez),3);
        
        for m = 1:length(kspace)
            neurordm        = squareform(dist(neurom(kspace(m),f)));
            subrdm          = squareform(eegrdm(:,:,s));
            corrmap(m,s,c)  = rankCorr_Kendall_taua(neurordm',subrdm');
        end
        
        [Y(s,c) I(s,c)] = max(corrmap(:,s,c));
        allI = find(corrmap(:,s,c) == Y(s,c));
        
        if allI(1) < find(kspace == 1)
            I(s,c) = allI(end);
        end
    end
end

bestk.number = kspace(I(:,1))';
bestk.donkey = kspace(I(:,2))';

%% k correlations

if do.signif
    corrdat 	= log([bestk.behav, bestk.number, bestk.donkey]);
    [r,p]       = corr(corrdat,'Type','Spearman')
end

%%

clear modresp

% Load paths
Numbers_load;

normnum     = (stim.samples - min(stim.samples(:)))/(max(stim.samples(:)) - min(stim.samples(:))); % positive input values
catdat    	= sign(stim.category - 1.5); % 1 = orange, 2 = blue
X           = [normnum catdat];
Y         	= data.chosenCat - 1;
feat_space  = unique(normnum)';
gnorm       = 0;

% Get subject responses
for s = 1:params.nsubj
    
    fprintf('\nBest fit sub %d\n',params.submat(s));
    
    idx         = data.sub == params.submat(s) & Y >= 0;
    frame       = unique(data.frame(idx));
    
    % Get response probabilities per number
    tmpSubResp  = makeLong(repmat(data.chosenCat(idx)-1,1,params.nsamp)');
    tmpCat      = makeLong(stim.category(idx,:)');
    tmpSamp     = makeLong(stim.samples(idx,:)');
    
    for c = 1:params.nstim
        subresp(s,c) = mean([tmpSubResp(tmpSamp == c & tmpCat == 2); -(tmpSubResp(tmpSamp == c & tmpCat == 1))+1]);
    end
    
end

% Subject performance
subresp_flip = subresp;
subresp_flip(params.frame == 0,:) = 1-subresp_flip(params.frame == 0,:);

% Get model responses
[~, pred] = psymodfun(median(bestfit),Y,X,1,params.nsamp,feat_space,gnorm,1);

% Get response probabilities per number
tmpModResp  = makeLong(repmat(pred,1,params.nsamp)');
tmpCat      = makeLong(stim.category');
tmpSamp     = makeLong(stim.samples');

for c = 1:params.nstim
    modresp(c) = mean([tmpModResp(tmpSamp == c & tmpCat == 2); -(tmpModResp(tmpSamp == c & tmpCat == 1))+1]);
end

% Fit k distortion to scale of behavioural data
nummod = [0:.2:1].^median(bestk.number);
nummod = (max(median(subresp_flip),[],2) - min(median(subresp_flip),[],2)) .* nummod + min(median(subresp_flip),[],2);

donkmod = [0:.2:1].^median(bestk.donkey);
donkmod = (max(median(subresp_flip),[],2) - min(median(subresp_flip),[],2)) .* donkmod + min(median(subresp_flip),[],2);

%% Plot figures

histcol     = colz([3,1,2],:);
xlims       = [-2 2];
nmax        = [0 1];
nbins       = 5;
titlez      = {'Behavior','EEG Numerical','EEG Bandit'};

figH = figure;
set(figH, 'Position',  [50, 50, 900, 400])

plotdat 	= log([bestk.behav, bestk.number, bestk.donkey]);

for h = 1:3
    subplot(2,3,h); hold on;
    fd      = fitdist(plotdat(:,h),'Kernel','BandWidth',.25);
    xvals   = linspace(-2.2,2.2,100);
    pdfY    = pdf(fd,xvals);
    hA      = area(xvals,pdfY,'FaceColor',histcol(h,:),'EdgeColor',histcol(h,:),'LineWidth',3,'FaceAlpha',.3);
    plot([0 0],nmax,'k--','LineWidth',2);
    plot([-2.5 2.5],[0 0],'Color',[1 1 1],'LineWidth',2.8);
    
    % Plot human distribution
    sk      = round(plotdat(:,h),1);
    kvals   = unique(sk);
    for s = 1:length(kvals)
        nk = length(find(sk == kvals(s)));
        plot([kvals(s) kvals(s)],[0 .025*nk],'k','LineWidth',1.5);
    end
    
    ax = gca;
    set(ax ,'Layer', 'Top');
    xlim(xlims);
    ylim(nmax);
    set(ax,'FontSize',axfntsz);
    xlabel('Log(k)','FontSize',axfntsz);
    if h == 1
        ylabel('P(k)','FontSize',16);
    end
    
    title(titlez{h},'FontWeight','bold');
    
end

for l = 1:3
    subplot(2,3,l+3); hold on;
    
    if l == 1
        plotdat = modresp;
    elseif l == 2
        plotdat = nummod;
    elseif l == 3
        plotdat = donkmod;
    end
    
    plot([0 params.nstim+1],[.5 .5],'k:','LineWidth',1.5);
    h(1) = plot(1:params.nstim,plotdat,'-','Color',histcol(l,:),'LineWidth',3);
    h(2) = plot(1:params.nstim,median(subresp_flip),'o','MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',2.5,'MarkerSize',8);
    xlim([.5 params.nstim+.5]);
    ylim([.35 .7]);
    ax = gca;
    set(ax,'XTick',1:params.nstim,'XTickLabel',{1:6},'FontSize',axfntsz);
    ax.XAxis.FontWeight = 'bold';
    if l == 1
        ylabel('P(choice)','FontSize',16);
    end
    
    legend(h,{'Model','Human'},'Location','NorthWest');
    legend boxoff
    
    % Change position
    ax.Position(1) = ax.Position(1)+.025;
    ax.Position(3) = ax.Position(3)-.05;
    
end

