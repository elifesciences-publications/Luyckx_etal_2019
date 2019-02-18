function [rdmSet] = rdm_computeRDM_crossval(results,disttype)
%% RDM_COMPUTEMAHALRDM
%
% computes distance RDMset
%
% results.betas  = electrode-time-condition matrix
% results.resids = electrode-time-trial matrix [for Mahalanobis]
% disttype = 'euclidean' / 'pearson' / 'spearman'

% we want a time-condition-condition rdm-set:
rdmSet = zeros(size(results.betas,3),size(results.betas,3),size(results.betas,2),size(results.betas,2));

% iterate through all time points
for t1 = 1:size(results.betas,2)
    for t2 = 1:size(results.betas,2)
        
        respMat1 = squeeze(mean(results.betas(:,t1,:),2));
        respMat2 = squeeze(mean(results.betas(:,t2,:),2));
        
        switch disttype
            case 'euclidean'
                rdmSet(:,:,t1,t2) 	= pdist2(respMat1',respMat2');
            case 'pearson'
                rdmSet(:,:,t1,t2)   = pdist2(respMat1',respMat2','correlation');
        end
    end
end

end