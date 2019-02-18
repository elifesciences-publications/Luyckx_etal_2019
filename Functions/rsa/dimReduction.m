function [xdat_reduced,Smat] = dimReduction(xdat,removeDim)
%function xdat_reduced = dimReduction(xdat,removeDim)
% 
% xdat: conditions x channels x timepoints
% removeDim: which dimension(s) to remove

xdat_reduced = 0*xdat;
Smat         = zeros(size(xdat,1),size(xdat,1),size(xdat,3));

for t = 1:size(xdat,3)
    
    respMat     = xdat(:,:,t);    
    [U,S,V]     = svd(respMat);
    Smat(:,:,t) = S(:,1:size(S,1));
        
    for i = removeDim
        S(i,i) = 0;
    end
    
    xdat_reduced(:,:,t) = U*S*V';
end

end

