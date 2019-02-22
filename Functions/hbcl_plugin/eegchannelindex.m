function [chann, EEG] = eegchannelindex( EEG, tempchann)
%   Looks up the channel index for the requested channel. Outputs
%   the channel index value or 0 if channel does not exist.
%
%   1   Input EEG File From EEGLAB
%   2   Channel label
%
%   index = eegchannelindex( EEG, 'FZ');
%
%   Author: Matthew B. Pontifex, Health Behaviors and Cognition Laboratory, Michigan State University, March 20, 2014
    chann = 0;
    n = size(EEG.chanlocs, 2);
    for m=1:n
        tempval = EEG.chanlocs(m).('labels');
        if (strcmp(tempval,tempchann) > 0)
            chann = m;
            break;
        end
    end
end
