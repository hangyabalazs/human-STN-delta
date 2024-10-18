function [EEG Lnfilt] = rem_line_noise(EEG)
% Removes power line noise (50 Hz) from EEG data (eeglab structure)
% First tries CleanLine, if noise has not been removed sufficiently (has to
% be checked visually on the appeared PSD), notch filter is applied (45-55
% Hz notch filter). Label of applied filter is stored in Lnfilt variable.

% Is there any line noise?
[fig1, fig2] = check_eegdata(EEG);
inp0 = input('Linenoise? If no linenoise, press 0, otherwise any key.\n');
% inp0 = 1;
%
close(fig1); close(fig2);


if inp0==0
    Lnfilt = 'NoLN';
    
else
    % CLEANLINE FILTER
    
    EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:EEG.nbchan] ,...
        'computepower',1,'linefreqs',50,'newversion',0,...
        'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,...
        'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',4);
    Lnfilt = 'CLN';
    
    [fig1, fig2] = check_eegdata(EEG, strcat('CleanLine filtered data'));
    
    inp = input('Linenoise removed? If yes, press 1, otherwise any key\n');
    %     inp = 0;
    
    if inp~=1
        close(fig1,fig2)
        % notch filter to remove 50 Hz linenoise
        EEG = pop_eegfiltnew(EEG, 'locutoff',45,'hicutoff',55,'revfilt',1,'plotfreqz',1);
        Lnfilt = 'notch';
        [fig1, fig2] = check_eegdata(EEG,'');
    end
    
    close(fig1,fig2);
    
    
end

%--------------------------------------------------------------------------
function [fig1, fig2] = check_eegdata(EEG, figtitle)
% Plots power spectrum

narginchk(1,2)
if nargin<2
    figtitle = '';
end

% plot eeg time series

pop_eegplot(EEG,1,1,1,[],'srate',EEG.srate,'spacing',75,...
    'eloc_file',EEG.chanlocs, 'winlength',5,'dispchans',32,'events',EEG.event,...
    'plottitle', figtitle);
fig1 = gcf;

% plot power spectrum

fig2 = figure;
pop_spectopo(EEG,1,[1,size(EEG.data,2)],'EEG','freqrange',[1 100],'title',figtitle);