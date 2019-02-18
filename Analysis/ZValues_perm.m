function [Z] = ZValues_perm(s,do,params,paths)
%function ZValues_perm(s,do,params,paths)

makeLong = @(x) x(:);

% Obtain model RDMs
mods    = ModelRDM;
models  = fieldnames(mods);

subvec  = setdiff(1:params.nsubj,s);
submat2 = params.submat(subvec);

%% Crossvalidation within participant

fprintf('\nCorrelating model - EEG RDM subject %d, %s phase.\n',params.submat(s),params.whichPhase);

% Load data
inputfile = sprintf('Crossval_%03d_RDM_%s_%s',params.submat(s),params.whichPhase,params.disttype);
load(fullfile(paths.data.save,inputfile));

modcont = zeros(length(rdm.timepoints),length(rdm.timepoints));

% Smoothing RDM
if do.smooth
    wdwsz    = 60/4; % size convolution kernel
    rdm.data = smoothRDM(rdm.data,wdwsz);
end

% EEG - Model correlation
actmod  = makeLong(mods.(models{1}));
for t1 = 1:length(rdm.timepoints)
    for t2 = 1:length(rdm.timepoints)
        rdmvec = makeLong(rdm.data(7:12,1:6,t1,t2));
        modcont(t1,t2) = rankCorr_Kendall_taua(actmod,rdmvec);
    end
end

% Baseline correct
tmp         = modcont;
tmp(rdm.timepoints >= 0, rdm.timepoints >= 0) = nan;
baseav      = nanmean(tmp(:));
modcont     = modcont - repmat(baseav,length(rdm.timepoints),length(rdm.timepoints));

%% Crossvalidation between participants

% Subject correlations
for z = 1:params.nsubj-1
    
    fprintf('\nCorrelating model - EEG RDM subject %d to subject %d.\n',params.submat(s),submat2(z));
    
    % Load data
    inputfile = sprintf('Crossval_%03d_RDM_%s_%s_perm%d',params.submat(s),params.whichPhase,params.disttype,submat2(z));
    load(fullfile(paths.data.save,inputfile));
    
    if z == 1
        permdat = zeros(length(submat2),length(rdm.timepoints),length(rdm.timepoints)); % initialise permutated data
    end
   
    % Smoothing RDM
    if do.smooth
        wdwsz    = 60/4; % size convolution kernel
        rdm.data = smoothRDM(rdm.data,wdwsz);
    end
    
    % EEG - Model correlation
    actmod  = makeLong(mods.(models{1}));
    for t1 = 1:length(rdm.timepoints)
        for t2 = 1:length(rdm.timepoints)
            rdmvec = makeLong(rdm.data(7:12,1:6,t1,t2));
            permdat(z,t1,t2) = rankCorr_Kendall_taua(actmod,rdmvec);
        end
    end
end

% Baseline correct
tmp         = permdat;
tmp(:,rdm.timepoints >= 0, rdm.timepoints >= 0) = nan;
baseav      = nanmean(reshape(tmp,[length(submat2),length(rdm.timepoints)^2]),2);
permdat     = permdat - repmat(baseav,1,204,204);

%% Calculate z-scores per subject

Z = (modcont - squeeze(nanmean(permdat)))./(squeeze(nanstd(permdat,[],1))./sqrt(params.nsubj-1));

% Save data
if do.save_zval
    save(fullfile(paths.data.save,sprintf('Zval_%03d_perm',params.submat(s))),'Z');
    fprintf('\nData subject %d saved.\n',s);
end

end

