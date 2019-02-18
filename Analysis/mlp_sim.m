function [loss,in1,in2,Whid] = mlp_sim(t,params,paths)
%function [loss,in1,in2,Whid] = mlp_sim(t,params,paths)

%% Initialise

loss    = zeros(2,2,params.iterations*2);
aloss   = loss;
Whid    = zeros(2,2,params.ninputs*2+1,params.nhid);
in1     = zeros(2,2,params.ninputs,params.nstim);
in2     = in1;

%%

fprintf('\nITERATION %d/%d\n',t,params.nsim);

% define random inputs
rr2= makeflip2(params.ninputs,params.noutputs,params.nstim);

% define random uniform supervision signal
tr2 = rand(params.noutputs,params.nstim)';

for expe = 1:2
        
    if expe==1  % A = shuffled
        
        %fprintf('\nA = random\n');
        
        [rr1, tr1] = makeflip2(params.ninputs,params.noutputs,params.nstim);
        rr1 = reshape(Shuffle(rr1(:)),size(rr1));
        tr1 = Shuffle(tr2);
        
    else  % inputA = inputB
        
        %fprintf('\nA = B\n');
        
        [rr1, tr1] = makeflip2(params.ninputs,params.noutputs,params.nstim);
        tr1 = tr2;
    end
    
    for shuf = 1:2  % shuffled weights, then unshuffled
        
%         if shuf == 1
%             fprintf('\nHidden weights shuffled\n');
%         else
%             fprintf('\nHidden weights unshuffled\n');
%         end
        
        % XOR input task 1
        XX1 = [ones(1,params.nstim);rr1;zeros(params.ninputs,params.nstim)]';  % inputs, including bias added manually
        
        % XOR input task 2
        XX2 = [ones(1,params.nstim);zeros(params.ninputs,params.nstim);rr2]';  % inputs, including bias added manually
        
        % Initialize the weights from input to hidden
        W1 = randn(size(XX1,2),params.nhid)*0.01;
        
        % initialize the weights from hidden to output (+1 for bias)
        W2 = randn(params.nhid,params.noutputs)*0.01;
        
        % run network on task 1, evaluate on task 2
        %fprintf('... training');
        [loss1,W1,W2,aloss1] = mlp(XX1,tr1,W1,W2,params.alpha,params.iterations,XX2,tr2);
        
        % shuffle weights for layer 2
        if shuf == 1
            W2 = reshape(Shuffle(W2(:)),size(W2,1),size(W2,2));
        end
        
        % run network again on task 2, evaluate on task 1
        %fprintf('... testing');
        [loss2,W1,W2,aloss2] = mlp(XX2,tr2,W1,W2,params.alpha,params.iterations,XX1,tr1);
        
        % save loss
        loss(expe,shuf,:)         = [loss1 loss2];
        aloss(expe,shuf,:)        = [aloss1 aloss2];
        in1(expe,shuf,:,:)        = rr1;
        in2(expe,shuf,:,:)        = rr2;
        Whid(expe,shuf,:,:)       = W1;
    end
end

%% Save iteration

save(fullfile(paths.data.save,sprintf('mlp_sim_run%d',t)),'loss','in1','in2','Whid');

end

