function [loss,W1,W2,loss2,sumW1,sumW2,sumd1,sumd2] = mlp(X,T,W1,W2,alpha,iterations,X2,T2)
% function [loss,W1,W2,loss2,sumW1,sumW2] = mlp(X,T,W1,W2,alpha,iterations,X2,T2);
%
% generic 3-layer MLP function
%
% inputs:
% X - input data [nstim x ninputs]
% T - supervision signal [nstim x noutputs]
% W1 - weight initialisation for layer1 -> 2 [ninputs x nhid]
% W1 - weight initialisation for layer2 -> 3 [nhid x noutputs]
% alpha - learning rate
% iterations - n training cycles
% X2 - auxiliary input, evaluated but not trained
% T2 - auxiliary supervision signal
%
% outputs:
% loss - loss
% W1 - first level weights, trained
% W2 - second level weights, trained
% loss2 - auxiliar loss
% sumW1 - Frobenius norm for trained weights level 1, per cycle
% sumW2 - Frobenius norm for trained weights level 2, per cycle

%% Functions
squelch = @(x) x(:);

%%

for i = 1:iterations
    
    sumW1(i) = sqrt(sum(abs(W1(:))));
    sumW2(i) = sqrt(sum(abs(W2(:))));
    
    % first evaluate on auxiliary data (e.g. test)
    
    % input to hidden
    H   = X2*W1;                       % hidden layer activations
    fH  = 1./(1+exp(-H));            % sigmoid activation function
    
    % hidden to output
    Y   = fH*W2;                      % output layer activations
    fY  = 1./(1+exp(-Y));            % sigmoid activation function
    
    % compute error
    E           = sum((T2-fY).^2,2);           % euclidean error
    loss2(i)    = sum(E,1);             % save loss is squared error
    
    % now on initial data (e.g. training)
    
    % input to hidden
    H   = X*W1;                       % hidden layer activations
    fH  = 1./(1+exp(-H));            % sigmoid activation function
    
    % hidden to output
    Y   = fH*W2;                      % output layer activations
    fY  = 1./(1+exp(-Y));            % sigmoid activation function
    
    % compute error
    E       = sum((T-fY).^2,2);           % euclidean error
    loss(i) = sum(E,1);             % save loss is squared error
    
    % BACK-PROPAGATION
    deriv_loss_loss = 1;
    deriv_loss_E    = repmat(deriv_loss_loss,size(E));
    deriv_loss_fY   = repmat(deriv_loss_E,[1,size(fY,2)]) .* 2 .* (fY - T);
    deriv_loss_Y    = deriv_loss_fY .* (1 - fY) .* fY;
    deriv_loss_fH   = deriv_loss_Y * W2';
    deriv_loss_W2   = fH' * deriv_loss_Y;
    deriv_loss_H    = deriv_loss_fH .* (1 - fH) .* fH;
    deriv_loss_X    = deriv_loss_H * W1';
    deriv_loss_W1   = X' * deriv_loss_H;
    
    % update weights
    W1 = W1 - (alpha * deriv_loss_W1); % update weights
    W2 = W2 - (alpha * deriv_loss_W2); % update weights
    
    sumd1(i) = sum(squelch(abs(deriv_loss_W1)));
    sumd2(i) = sum(squelch(abs(deriv_loss_W2)));
    
end

