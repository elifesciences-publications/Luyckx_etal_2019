function eegplugin_hbcl(fig,try_strings,catch_strings)

version = 1.0;
% create menu
toolsmenu1 = findobj(fig, 'tag', 'tools');
toolsmenu2 = findobj(fig, 'tag', 'ERPLAB');

uimenu( toolsmenu1, 'label', 'Visually Inspect EEG - Plot Array', 'separator','on','callback', '[EEG LASTCOM]=pop_visualinspectarray(EEG);  [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);');
uimenu( toolsmenu1, 'label', 'Visually Inspect EEG - Plot Trials', 'callback', '[EEG LASTCOM]=pop_visualinspecttrials(EEG); [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);');

uimenu( toolsmenu2, 'label', 'Measure Peak Interval', 'separator','on','callback', '[ERP_MEASURES]=pop_erppeakinterval(ERP);');
uimenu( toolsmenu2, 'label', 'Egg Head Plot', 'callback', 'pop_eggheadplot');
uimenu( toolsmenu2, 'label', 'Plot ERP Array', 'callback', '[LASTCOM]=pop_plotERParray(ERP);');
uimenu( toolsmenu2, 'label', 'Compare ALLERP Arrays', 'callback', '[LASTCOM]=pop_compareERParray(ALLERP);');
