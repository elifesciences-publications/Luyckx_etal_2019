function [bestparam,mfit,moddata] = RLmodel(s,data,stim,use,init_param,init_prob,crossval)
%function [bestparam,mfit,moddata] = RLmodel(s,data,stim,use,init_param,init_prob,crossval)
%   Based on juechems_etal_2017 (Github)

nparams = size(init_param,1); % number of parameters
bnd_low = min(init_param,[],2)'; % lower bound
bnd_up  = max(init_param,[],2)'; % upper bound
submat  = unique(data.sub)'; % subject matrix
ntrials = sum(use); % number of valid trials

init    = rand(1,nparams);
init    = init .* (bnd_up - bnd_low) + bnd_low;

% Find global minimum
gs                  = GlobalSearch('Display','off');
problem             = createOptimProblem('fmincon','x0',init,'objective',@log_likelihood,'lb',bnd_low,'ub',bnd_up);
[bestparam,mfit]    = run(gs,problem);

%% CROSS-VALIDATION

if crossval
    use     = data.sub == submat(s) & ~use;
    ntrials = sum(use);
    mfit    = log_likelihood(bestparam);
end

% Get model choices and probabilities
moddata = RLmodel(bestparam);

%% LOG LIKELIHOOD FUNCTION

    function log_l = log_likelihood(parms)
        mdata           = RLmodel(parms);
        decision        = data.chosenOrder(use) - 1;
        valid_tr        = decision >= 0;
        
        mdata.p_choice(mdata.p_choice < 0.0000001) = 0.0000001;
        mdata.p_choice(mdata.p_choice > 0.9999999) = 0.9999999;
        
        likelihood      = decision(valid_tr) .* log(mdata.p_choice(valid_tr)) + (1-decision(valid_tr)) .* log(1-mdata.p_choice(valid_tr));
        log_l           = -1*sum(likelihood);
    end

%% RL MODEL

    function mdata = RLmodel(parms)
        
        mdata.p_ranks       = zeros(ntrials,length(unique(stim.combo)));
        mdata.p_choice      = nan(ntrials,1);
        
        combos              = stim.combo(use,:); % bandit options on trial t
        subresp             = data.chosenOrder(use,:); % ptt choice first or second?
        outcomes            = data.outcome(use); % feedback
        chosen_rank         = data.chosenProb(use); % which bandit chosen
        mdata.chosen_order  = 0*chosen_rank;
        mdata.chosen_bandit = mdata.chosen_order;
        mdata.p_ranks(1,:) 	= init_prob(s,:); % inital values bandits
        
        for t = 1:ntrials
            
            % Current estimated probabilities
            pnow 	= mdata.p_ranks(t,:);
            dv  	= pnow;
            
            % Model probability
            mdata.p_choice(t) = sigmoid(diff(dv(combos(t,:))),0.05,0.9,0,parms(2));
            
            if mdata.p_choice(t) == .5
                mdata.chosen_order(t) = randi(2);
            else
                mdata.chosen_order(t) = .5*sign(mdata.p_choice(t)-.5)+1.5;
            end
            
            mdata.chosen_bandit(t) = combos(t,mdata.chosen_order(t));
            
            % Update ranks according to participant's choice
            if t < ntrials
                mdata.p_ranks(t+1,:) = pnow;
                if subresp(t) > 0 % when ppt didn't give response, no update
                    mdata.p_ranks(t+1,chosen_rank(t)) = pnow(chosen_rank(t)) + parms(1)*(outcomes(t) - pnow(chosen_rank(t)));
                end
            end
            
        end
        
    end
end
