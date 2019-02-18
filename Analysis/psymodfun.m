function [G, pred, g] = psymodfun(b,Y,X,ML,nk,f,gnorm,frame)
%
% logfunX_ml_power
% psychometric model for sequential comparison tasks
%
% Input:
%  	- b: [motor bias,kappa,late noise,leak];
%	- Y: subject choices
% 	- X: [number values, category index (-1 1)]
% 	- ML: use maximum likelihood estimation (FALSE/[TRUE])
%	- nk: number of samples
% 	- f: feature space (-1 -> 1)
%  	- gnorm: normalise function by g (FALSE/[TRUE])
%   - frame: flips for low frame
%
%  Output:
%  	- G: fit of model (max likelihood)
%  	- pred: predicted choices
% 	- g: normalization factor,
%
% adapted from (c)Bernhard Spitzer, 2016

sz = length(b);

% weighting function
pownum = abs((X(:,1:nk)).^b(2)); %b(2): kappa

% g  normalization
g = sum(abs(f).^b(2)) / sum(abs(f));

if gnorm
    pownum=pownum/g;
end

% sign-flip according to sample category;
pflip   = pownum.*X(:,nk+1:2*nk);

% Leak
leaker  = (1-b(4)).^(nk-[1:nk]'); % leakage term
pfin    = pflip*leaker; % if b(4)=0, this corresponds to sum(pflip,2)

% Choice
DV      = b(1)+pfin; % decision value
pred    = sigmoid(DV,0,1,0,b(3)); %1./(1+exp(-DV./b(3))); % logit rule - b(3): sigma (noise)

if frame == 0
    pred = 1-pred;
end

% ML fit (not used in simulation)
if ML == 1
    G  	= -1*sum(Y.*log(pred) + (1-Y).*log(1-pred));
else
    G 	= sum((Y-pred).^2); % ols
end


end