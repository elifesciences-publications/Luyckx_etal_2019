function [rdmSet,whiteEEG] = rdm_computeRDM(results,disttype)
%% RDM_COMPUTEMAHALRDM
%
% computes distance RDMset
%
% results.betas  = electrode-time-condition matrix
% results.resids = electrode-time-trial matrix [for Mahalanobis]
% disttype = 'euclidean' / 'correlation' / 'mahalanobis'
%
% (c) Timo Flesch, 2016

% we want a time-condition-condition rdm-set:
rdmSet = zeros(size(results.betas,2),size(results.betas,3),size(results.betas,3));

% iterate through all time points
for timePoint = 1:size(results.betas,2)
    
    respMat     = squeeze(mean(results.betas(:,timePoint,:),2));
    
    switch disttype
        case 'euclidean'
            rdmSet(timePoint,:,:)   = squareform(pdist(respMat'));
        case 'mahalanobis'
            residMat                = squeeze(mean(results.resids(:,timePoint,:),2));
            rdmSet(timePoint,:,:)   = squareform(pdist(respMat','mahalanobis',covdiag(residMat')));
        case 'pearson'
            %residMat                = squeeze(mean(results.resids(:,timePoint,:),2));
            %whiteEEG(timePoint,:,:) = respMat'*covdiag(residMat')^(-.5);
            rdmSet(timePoint,:,:)   = squareform(pdist(respMat','correlation'));           
    end
        
end

end