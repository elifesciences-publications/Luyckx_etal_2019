function varargout = meeze(varargin)
    varargout{:} = squeeze(nanmean(varargin{:}));
end