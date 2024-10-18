function EEG_ep = substr_bas(EEG_ep,baseline_win)


% Subtract baseline from each channel and epoch

% EEG_ep: epoched EEG data (eeglab format)
% baseline_win: baseline window relative to event timetamps in seconds (ex: [-2 -1])
%-----------------------------------------------
bas_fr = abs(baseline_win(1,1)) + abs(baseline_win(1,2))*EEG_ep.srate;
EEG_ep.data = rmbase(EEG_ep.data,[],[1:bas_fr]); % dataset 1
