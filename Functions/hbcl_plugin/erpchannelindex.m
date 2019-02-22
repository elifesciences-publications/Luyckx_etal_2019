function [chann, ERP] = erpchannelindex( ERP, tempchann)
%   Looks up the index number for the specified channel
%
%   [ERP, chann] = erpchannelindex(ERP, 'M1');

    chann = 0;
    n = size(ERP.chanlocs, 2);
    for m=1:n
        tempval = ERP.chanlocs(m).('labels');
        if (strcmp(tempval,tempchann) > 0)
            chann = m;
            break;
        end
    end
end
