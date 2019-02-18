function mods = ModelRDM_CPP(s,params)
%(c) Bernhard Spitzer, 2016

%% Numerical distance

pc=zeros(6,6);
for i=1:6
    for j=1:6
        pc(i,j)=abs(i-j)./5;
    end
end

mods.numd=pc;

%% CPP Numerical

savefolder = ''; % unused
figfolder = ''; % unused
Numbers_load;
load(fullfile(paths.data.EEG.raw,'Peak_CPP_numbers_collapsed'));

bardat_num     	= peakz(s,:);
mods.numCPP     = dist(bardat_num);

%% CPP Bandit

savefolder = ''; % unused
figfolder = ''; % unused
Bandit_load;
load(fullfile(paths.data.EEG.raw,'Peak_CPP_donkey_collapsed'));

bardat_donkey   = peakz(s,:);
mods.donkeyCPP  = dist(bardat_donkey);

%% Cross-validation

crossCPP        = dist([bardat_num,bardat_donkey]);
mods.crossCPP   = crossCPP(7:12,1:6);

end