function smoothdat = smoothRDM(rdmdat,wdwsz)
% function smoothdat = smoothRDM(rdmdat,wdwsz)
%
% rdmdat    = condition x condition x time (x time)
% wdwsz     = n data points (not time points)
%
% smoothdat = smoothed RDM

fprintf('... smoothing\n');

datsize     = ndims(rdmdat);
smoothdat   = 0*rdmdat;

if datsize == 3    
    smoothfilt  = unifpdf(1:wdwsz,1,wdwsz);
    
    for i = 1:size(rdmdat,1)
        for j = 1:size(rdmdat,2)
            smoothdat(i,j,:) = conv(squeeze(rdmdat(i,j,:)),smoothfilt,'same');
        end
    end
           
elseif datsize == 4
    smoothfilt  = fspecial('average',wdwsz);
    smoothdat   = rdmdat*0;
    
    for i = 1:size(rdmdat,1)
        for j = 1:size(rdmdat,2)
            smoothdat(i,j,:,:) = conv2(squeeze(rdmdat(i,j,:,:)),smoothfilt,'same');
        end
    end
end

end

