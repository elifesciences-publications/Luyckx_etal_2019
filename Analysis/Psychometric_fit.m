%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 3A: psychometric model fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DATA

clc
clear

% Paths
savefolder  = 'Numbers_model'; % folder to save newly created data (unused)
figfolder   = 'Behaviour';

% Load stuff
Numbers_load; % load data and path
Config_plot; % load plot variables

% Logicals
do.modelfit     = true; % do model fitting?
do.save_fit     = true; % save model estimates?
do.plotting     = true; % plot?
do.signif       = false; % test significance?

%% Extra variables

normnum     = (stim.samples - min(stim.samples(:)))/(max(stim.samples(:)) - min(stim.samples(:))); % positive input values
catdat    	= sign(stim.category - 1.5); % 1 = orange, 2 = blue
X           = [normnum catdat];
Y         	= data.chosenCat - 1;
feat_space  = unique(normnum)';

% Initialise parameters
offStart    = 0; % bias b
kappaStart  = 1; % kappa
noiseStart  = 1; % late noise
leakStart   = 0; % leak

kmax        = 10;
smax        = 8;
letLeak     = 1;
b0          = [0 kappaStart noiseStart leakStart]; % initialise values
lb          = [-Inf 0.1 0.001 0]; % lower bound
ub          = [Inf kmax smax letLeak]; % upper bound
gnorm       = 0;

optimist    = optimset('MaxFunEvals',50000,'MaxIter',50000,'Display','off');

%% Run parameter estimation

if do.modelfit
    % Fit model
    for s = 1:params.nsubj
        
        fprintf('\nFitting sub %d\n',params.submat(s));
        
        idx     = data.sub == params.submat(s) & Y >= 0;
        frame   = unique(data.frame(idx));
        
        [bestfit(s,:),LL(s)] = fmincon(@(b) psymodfun(b,Y(idx),X(idx,:),1,params.nsamp,feat_space,gnorm,frame),b0,[],[],[],[],lb,ub,[],optimist);
    end
else
    % Load fitted data
    load(fullfile(paths.data.model,sprintf('Numbers_modelfit')));
end

% Get model responses
for s = 1:params.nsubj
    
    fprintf('\nBest fit sub %d\n',params.submat(s));
    
    idx         = data.sub == params.submat(s) & Y >= 0;
    frame       = unique(data.frame(idx));
    
    numObs(s)   = sum(data.sub == params.submat(s) & data.r > 0);
    
    [Gopt(s,:), pred, g(s)] = psymodfun(bestfit(s,:),Y(idx),X(idx,:),1,params.nsamp,feat_space,gnorm,frame);
    
    % Get response probabilities per number
    tmpModResp  = makeLong(repmat(pred,1,params.nsamp)');
    tmpSubResp  = makeLong(repmat(data.chosenCat(idx)-1,1,params.nsamp)');
    tmpCat      = makeLong(stim.category(idx,:)');
    tmpSamp     = makeLong(stim.samples(idx,:)');
    
    for c = 1:params.nstim
        modresp(s,c) = mean([tmpModResp(tmpSamp == c & tmpCat == 2); -(tmpModResp(tmpSamp == c & tmpCat == 1))+1]);
        subresp(s,c) = mean([tmpSubResp(tmpSamp == c & tmpCat == 2); -(tmpSubResp(tmpSamp == c & tmpCat == 1))+1]);
    end
    
end

if do.save_fit
    save(fullfile(paths.data.save,sprintf('Numbers_modelfit')),'bestfit','modresp');
    fprintf('\n %s saved.\n',sprintf('Numbers_modelfit'));
end

%% Plot
if do.plotting
    
    % Average performance vs model
    subresp_flip = subresp;
    subresp_flip(params.frame == 0,:) = 1-subresp_flip(params.frame == 0,:);
    modresp_flip = modresp;
    modresp_flip(params.frame == 0,:) = 1-modresp_flip(params.frame == 0,:);
    
    avsub = mean(subresp_flip);
    SEsub = std(subresp_flip)/sqrt(params.nsubj);
    avmod = mean(modresp_flip);
    SEmod = std(modresp_flip)/sqrt(params.nsubj);
    
    figB = figure; hold on;
    plot([0 params.nstim+1],[.5 .5],'k:','LineWidth',1.5);
    h(1) = myfillsteplot(1:params.nstim,modresp_flip,colz(3,:),2);
    h(2) = plot(1:params.nstim,avsub,'ko','MarkerFaceColor','k','LineWidth',2,'MarkerSize',8);
    
    for e = 1:params.nstim
        plot([e e],[avsub(e)-SEsub(e) avsub(e)+SEsub(e)],'k','LineWidth',2);
    end
    
    legend(h,{'Model','Human'},'FontSize',titlefntsz,'Location','NorthWest');
    legend boxoff
    ax = gca;
    set(ax,'XTick',1:params.nstim,'XTickLabel',{1:6},'FontSize',axfntsz);
    ax.XAxis.FontWeight = 'bold';
    ax.XAxis.FontSize   = labfntsz;
    xlim([.5 params.nstim+.5]);
    ylim([.35 .7]);
    xlabel('Digit','FontSize',labfntsz);
    ylabel('Choice probability','FontSize',labfntsz);
end

