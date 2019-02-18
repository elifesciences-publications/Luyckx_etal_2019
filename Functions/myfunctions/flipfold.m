function flipfold_x = flipfold(xdat,zerodiag)
%function flipfold_x = flipfold(xdat,zerodiag)
%
% Flip and fold square matrix around diagonal
%
% xdat: square matrix
% zerodiag: zero the diagonal? [true]/false
%
% Fabrice Luyckx, 11/1/2019

    flipX       = rot90(fliplr(xdat));
    flipfold_x  = mean(cat(ndims(xdat)+1,xdat,flipX),ndims(xdat)+1);
    
    if zerodiag
        sz      = size(flipfold_x);
        nloops  = prod(sz(3:end));
        tmpX    = reshape(flipfold_x,[sz(1),sz(1),nloops]);
        zeroX   = zeros(sz(1),sz(1),nloops);
        
        for n = 1:nloops
            tmpn = tmpX(:,:,n);
            tmpn(logical(eye(sz(1)))) = 0;
            zeroX(:,:,n) = tmpn;
        end
        
        if length(sz) > 3
            flipfold_x = reshape(zeroX,sz(3:end));
        else
            flipfold_x = zeroX;
        end
    end
end

