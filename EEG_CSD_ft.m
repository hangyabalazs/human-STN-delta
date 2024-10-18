function EEG_CSD_ft(s2a)

% Applies CSD tranformation on EEG data 
% Uses CSD toolbox by Jürgen Kayser, 2009

% s2a: data to analyse (see getdata2analyse)

for snr = 1:length(s2a)
    
    try
    EEG = pop_loadset(fullfile(s2a(snr).folder,'EEG_2plot.set'));
    catch
        fprintf('No EEG\n');
        continue
%         pause
    end
    
    EEG = pop_currentdensity(EEG, 'method','spline');
    
    try
    pop_saveset(EEG,fullfile(s2a(snr).folder,'EEG_2plot_CSD.set'))
    catch
        pause
    end
end