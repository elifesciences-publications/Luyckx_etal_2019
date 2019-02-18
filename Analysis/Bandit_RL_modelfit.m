%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RL model to estimate bandit ranks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DATA

clc
clear

% Paths
savefolder          = 'Bandit_model'; % folder to save newly created data
figfolder           = 'Model'; % folder to save figures to
params.whichPhase   = 'test'; % use test phase data

% Load stuff
Bandit_load; % load data and path
Config_plot; % load plot variables

% Logicals
do.save_fit         = true; % save model fit?

%% Initialise values based on learning phase

% Load learning phase
params.whichPhase2 = 'learn';

% Load data
inputfile = sprintf('Donkey_fulldata_%s_all',params.whichPhase2);
load(fullfile(paths.data.behav,inputfile));

% Get chosen probability
data.chosenProb = zeros(length(data.sub),1);
init_prob    = zeros(params.nsubj,params.nsamp);

for t = 1:length(data.sub)
    if data.r(t) < 0
        data.chosenProb(t) = nan;
    else
        if data.r(t) == data.rmap(t)
            data.chosenProb(t) = stim.combo(t,1);
        else
            data.chosenProb(t) = stim.combo(t,2);
        end
    end
end

% Choice probabilities per bandit
for s = 1:params.nsubj
    idx = data.sub == params.submat(s) & data.r >= 0 & data.forcedChoice == 0;
    for i = 1:params.nstim
        tot = length(find(stim.combo(idx,:) == i));
        init_prob(s,i)   = length(find(data.chosenProb(idx) == i))/tot;
    end
end

%% Create necessary variables

% Load original phase again
inputfile = sprintf('Donkey_fulldata_%s_all',params.whichPhase);
load(fullfile(paths.data.behav,inputfile));

% Get chosen bandit
data.chosenProb     = zeros(params.ttrials,1);
data.chosenOrder    = data.chosenProb;

for t = 1:params.ttrials
    if data.r(t) < 0
        data.chosenProb(t)  = nan;
        data.chosenOrder(t) = nan;
    else
        if data.r(t) == data.rmap(t)
            data.chosenProb(t)  = stim.combo(t,1);
            data.chosenOrder(t) = 1;
        else
            data.chosenProb(t)  = stim.combo(t,2);
            data.chosenOrder(t) = 2;
        end
    end
end

%% Model fit

% Parameter bounds
init_param(1,:) = [.0001 .5]; % learning rate
init_param(2,:) = [.0001 .5]; % temperature

mod.bestparam       = [];
mod.LL              = [];
mod.param           = [];
mod.LL              = [];
mod.p_choice        = [];
mod.p_ranks         = [];
mod.ranks           = [];
mod.chosen_order    = [];
mod.chosen_bandit   = [];

for s = 1:params.nsubj
    fprintf('Fitting RL model, subject %d\n',params.submat(s));
    
    use = data.sub == params.submat(s);
    [mod.bestparam(:,s), mod.LL(s), moddata] = RLmodel(s,data,stim,use,init_param,init_prob,do.crossval);
    
    % Get rank indices
    [tmp, I] = sort(moddata.p_ranks,2);
    moddata.ranks = I;

    mod.p_choice        = [mod.p_choice; moddata.p_choice];
    mod.p_ranks         = [mod.p_ranks; moddata.p_ranks];
    mod.ranks           = [mod.ranks; moddata.ranks];
    mod.chosen_order    = [mod.chosen_order; moddata.chosen_order];
    mod.chosen_bandit   = [mod.chosen_bandit; moddata.chosen_bandit];
  
end

% Save concatenated data
if do.save_fit
    save(fullfile(paths.data.model,sprintf('Modelfit_full_%s_model_%s_global',params.whichPhase,curr_mod)),'mod');
    fprintf('\nData model %s saved.\n',curr_mod);
end
