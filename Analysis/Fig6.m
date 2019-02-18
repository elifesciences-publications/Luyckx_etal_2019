%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 6: neural network simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DATA

clc
clear

% Paths
paths.main              = fullfile(''); % change path to location of folder

if isempty(paths.main)
    error('Please change paths.main to the location of the folder ''Luyckx_etal_2019''.');
end

paths.simulation        = fullfile(paths.main,'Analysis');
paths.functions.main    = fullfile(paths.main,'Functions');
paths.data.main         = fullfile(paths.main,'Data','Bandit');
paths.data.save         = fullfile(paths.data.main,'NN'); % folder to save newly created simulation data
paths.figures.current   = fullfile(paths.main,'Figures','NN'); % folder to save plots

cd(paths.simulation);
addpath(paths.data.save);
addpath(genpath(paths.functions.main));

% Load stuff
Config_plot; % load plot variables

% Logicals
do.sim          = false;
do.rsa          = true;
do.signif       = true;
do.plotting     = true;

%% Neural network

% Network size
params.nstim        = 6; % number of stimuli per task
params.ninputs      = 20;% number of input units per module
params.noutputs     = 10; % number of output units
params.nhid         = 10; % number of hidden units

% Learning rate
params.alpha        = 0.001;

% Iterations of process
params.nsim         = 100;

% Number of learning iterations
params.iterations   = 1e6;

if do.sim
    
    for t = 1:params.nsim  % do this multiple times to get stable estimates
        mlp_sim(t,params,paths);
    end
    
else
    
    % Concatenate data    
    tot.loss    = zeros(params.nsim,2,2,params.iterations*2);
    tot.Whid    = zeros(params.nsim,2,2,params.ninputs*2+1,params.nhid);
    tot.in1     = zeros(params.nsim,2,2,params.ninputs,params.nstim);
    tot.in2     = tot.in1;
    
    for t = 1:params.nsim
        fprintf('\nLoading run %d ...\n',t);
        
        load(fullfile(paths.data.save,sprintf('mlp_sim_run%d',t)));
        
        tot.loss(t,:,:,:)   = loss;
        tot.Whid(t,:,:,:,:) = Whid;
        tot.in1(t,:,:,:,:)  = in1;
        tot.in2(t,:,:,:,:)  = in2;       
    end
    
    params.iterations = size(tot.loss,4)/2;
    
end

%% RSA

if do.rsa
    
    % Initialise
    hid_RSA1    = zeros(params.nsim,2,2,6,6);
    hid_RSA2    = hid_RSA1;
    hid_RSA12   = hid_RSA1;
    RDMcorr     = zeros(params.nsim,3,2,2);
    
    for t = 1:params.nsim
        
        fprintf('\nRSA run %d\n',t);
        
        for expe  = 1:2
            for shuf = 1:2
                for n = 1:params.nstim
                    XXt = [ones(1,1);squeeze(tot.in1(t,expe,shuf,:,n));zeros(params.ninputs,1)]';  % inputs, including bias added manually
                    hid_act1(n,:) = XXt*squeeze(tot.Whid(t,expe,shuf,:,:));                 % hidden layer activations
                    
                    XXt = [ones(1,1);zeros(params.ninputs,1);squeeze(tot.in2(t,expe,shuf,:,n))]';  % inputs, including bias added manually
                    hid_act2(n,:) = XXt*squeeze(tot.Whid(t,expe,shuf,:,:));                 % hidden layer activations
                end
                
                RSA1    = dist(hid_act1');
                RSA2    = dist(hid_act2');
                RSA12   = pdist2(hid_act1,hid_act2);
                
                hid_RSA1(t,expe,shuf,:,:)   = RSA1;
                hid_RSA2(t,expe,shuf,:,:)   = RSA2;
                hid_RSA12(t,expe,shuf,:,:)  = RSA12;
                
                % model RDM correlations
                modelRDM = squareform(pdist((1:6)'));
                
                RDMcorr(t,1,expe,shuf) = rankCorr_Kendall_taua(squareform(modelRDM),squareform(RSA1));
                RDMcorr(t,2,expe,shuf) = rankCorr_Kendall_taua(squareform(modelRDM),squareform(RSA2));
                RDMcorr(t,3,expe,shuf) = rankCorr_Kendall_taua(modelRDM(:),RSA12(:));
            end
        end
    end
end

%% Plots

if do.plotting
      
    %% Fig 6A: plot inputs

    nstim           = 6;
    nfeat           = 9;
    savenamez       = {'Xa','Xb'};
    
    for r = 1:2
        rr          = makeflip2(params.ninputs,params.noutputs,params.nstim);
        plotdat     = rr(1:nfeat,1:nstim);
        
        figI = figure;
        colormap(bone);
        h = imagesc(plotdat); hold on;
        
        ax = gca;
        set(ax,'xtick', linspace(0.5,nstim+0.5,nstim+1), 'ytick', linspace(0.5,nfeat+.5,nfeat+1));
        set(ax,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k');
        ax.GridAlpha = 1;
        
        set(ax,'XTick',ax.XTick,'XTickLabel',{''},'YTick',ax.YTick,'YTickLabel',{''});
        xlabel('STIMULUS','FontSize',labfntsz*3);
        ylabel('FEATURE','FontSize',labfntsz*3);
        
        if r == 2
            ax.YAxisLocation = 'right';
        end
        
        % RDM
        figA = figure;
        colormap(hot);
        imagesc(dist(rr));
        set(gca,'XTick',[],'YTick',[]);
        axis square;
        
        xlabel('STIMULUS','FontSize',labfntsz*3);
        ylabel('STIMULUS','FontSize',labfntsz*3);
        
        if r == 1
            ax.YAxisLocation = 'right';
        end

    end
    
    %% Fig. 6B: Loss over time

    titlez      = {'Shuffled','Unshuffled'};
    condlabz    = {'x_a = random','x_a , x_b: same structure'};
    firstits    = .2*params.iterations; % first iterations to plot 
    stepz       = .01*firstits; % downsample
    gapsz       = 15;
    rest        = 10;
    maxY        = 3;
       
    for h = 2%1:2
        figS(h) = figure; hold on;
                
        % Get subset data to plot
        traindat    = squeeze(tot.loss(:,:,h,1:firstits));
        plotdat     = traindat(:,:,1:stepz:end);
        cut         = nan(size(plotdat,1),size(plotdat,2),gapsz);
        plotdat     = cat(3,plotdat,cut);
        retraindat  = squeeze(tot.loss(:,:,h,params.iterations+[1:firstits]-stepz*rest));
        plotdat     = cat(3,plotdat,retraindat(:,:,1:stepz:end));
        plotdat     = cat(3,plotdat,cut);
        plotdat     = permute(plotdat,[1,3,2]);
        
        plot([size(plotdat,2)/2+rest size(plotdat,2)/2+rest],[0 maxY],'k:','LineWidth',2);
        
        steh(1) = plot(1:size(plotdat,2),squeeze(mean(plotdat(:,:,1))),'Color',nncolz(1,:),'LineWidth',4);
        steh(2) = plot(1:size(plotdat,2),squeeze(mean(plotdat(:,:,2))),'Color',nncolz(2,:),'LineWidth',4);

        % Axes
        ax = gca;
        set(ax,'FontSize',axfntsz);
        xlim([0 size(plotdat,2)]);
        ylim([0 maxY]);
        xlabel('Cycles (10^3)','FontSize',labfntsz);
        ylabel('Loss','FontSize',labfntsz);
        box off
        
        % Tick labels  
        ax.XTick    = [[1 50 100] [1 50 100]+gapsz+rest+100];
        xticklabz   = [1 100 200 1 100 1000];
        set(ax,'XTick',ax.XTick,'XTickLabel',xticklabz);
        
        if h == 2
            
            % Custom legend
            coord       = [160, 2.7];
            stepsz      = .4;
            for c = 1:length(condlabz)
                text(coord(1),coord(2)-(stepsz*(c-1)),condlabz{c},'Color',nncolz(c,:),'FontSize',28,'FontWeight','bold');
            end
            
        end  
    end 
    
    %% Fig 6C: retraining of shuffled
    
    firstits    = .2*params.iterations; % first iterations to plot
    stepz       = .01*firstits; % downsample
    gapsz       = 15;
    rest        = 10;
    maxY        = 3;
    
    figR = figure;
    hold on;
    
    % Get subset data to plot
    retraindat  = squeeze(tot.loss(:,:,1,params.iterations+[1:firstits]-stepz*rest));
    plotdat     = retraindat(:,:,1:stepz:end);
    cut         = nan(size(plotdat,1),size(plotdat,2),gapsz);
    plotdat     = cat(3,plotdat,cut);
    plotdat     = permute(plotdat,[1,3,2]);
    
    steh(1) = plot(1:size(plotdat,2),squeeze(mean(plotdat(:,:,1))),'Color',colz(1,:),'LineWidth',4);
    steh(2) = plot(1:size(plotdat,2),squeeze(mean(plotdat(:,:,2))),'Color',colz(2,:),'LineWidth',4);
    
    plot([rest rest],[0 maxY],'k:','LineWidth',2);
    
    % Axes
    ax = gca;
    set(ax,'FontSize',axfntsz);
    xlim([0 size(plotdat,2)]);
    ylim([0 maxY]);
    xlabel('Cycles (10^3)','FontSize',labfntsz);
    ylabel('Loss','FontSize',labfntsz);
    box off
    
    % Tick labels
    ax.XTick    = [0 rest+1 size(plotdat,2)/2 size(plotdat,2)];
    xticklabz   = {'...','1','100','1000'};
    set(ax,'XTick',ax.XTick,'XTickLabel',xticklabz);
    
    %% Fig 6D: summarise early loss
    
    indx        = 1+(size(tot.loss,4)/2):1e3+(size(tot.loss,4)/2);
    tmp         = mean(tot.loss(:,:,:,indx),4);
    testdat     = fliplr(reshape(tmp,[100,4])); % shuff,Ar/unshuff,Ar/shuff,AB/unshuff,AB
    plotdat     = mean(testdat);
    semdat      = std(testdat)/sqrt(params.nsim);
    testlabz    = fliplr({'Sh','Sh','Unsh','Unsh'});
    ncondz      = length(testlabz);
    barcolz     = fliplr([1 2 1 2]);
    
    % Repeated-measures ANOVA
    if do.signif
        t = table([1:params.nsim]',testdat(:,1),testdat(:,2),testdat(:,3),testdat(:,4));
        rm = fitrm(t,'Var2-Var5~Var1');
        ranova(rm)
    end
    
    % Plot
    figB = figure;  hold on;
    
    for b = 1:ncondz
        hB(b) = bar(b,plotdat(b),'FaceColor',nncolz(barcolz(b),:),'FaceAlpha',.8,'EdgeColor',nncolzedge(barcolz(b),:),'LineWidth',2);
    end
    
    errorbar(1:ncondz,plotdat,semdat,'k.','LineWidth',2);
    
    ax = gca;
    set(ax,'FontSize',axfntsz);
    set(ax,'XTick',[1.5 3.5],'XTickLabel',{'W2_{unshuff}','W2_{shuff}'});
    ax.XAxis.FontSize = labfntsz;
    ylabel('Loss','FontSize',labfntsz);
    box off
    xlim([.5 4.5]);
    ylim([0 5.1]);
    
    %% Fig 6E: RSA
    
    if do.rsa
        
        % Only A-B
        abdat       = fliplr(reshape(squeeze(RDMcorr(:,3,:)),[params.nsim 4]));
        plotdat     = mean(abdat);
        stedat      = std(abdat)./sqrt(params.nsim);
        testlabz    = fliplr({'Sh','Sh','Unsh','Unsh'});
        barcolz     = fliplr([1 2 1 2]);
        
        figAB = figure;  hold on;
        
        for b = 1:4
            hB(b) = bar(b,plotdat(b),'FaceColor',nncolz(barcolz(b),:),'FaceAlpha',.8,'EdgeColor',nncolzedge(barcolz(b),:),'LineWidth',2);
        end
        
        errorbar(1:4,plotdat,stedat,'k.','LineWidth',2);

        ax = gca;
        set(ax,'XTick',[1.5 3.5],'XTickLabel',{'W2_{unshuff}','W2_{shuff}'});
        ax.XAxis.FontSize = labfntsz;
        ylabel('Mean Kendall \tau_A','FontSize',labfntsz);
        box off
        xlim([.5 4.5]);
        ylim([0 .6]);
        
        % Test significance
        if do.signif
            [h,p,ci,stats] = ttest(abdat)
        end
    
        %% RSA crossvalidation rdm
        
        plotdat = squeeze(mean(hid_RSA12(:,2,2,:,:)));
        
        figR = figure;
        colormap(hot);
        imagesc(plotdat);
        set(gca,'XTick',[],'YTick',[]);
        axis square
        
    end
end
