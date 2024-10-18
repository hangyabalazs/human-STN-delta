function [EEG_ep1,EEG_ep2] =  prep_epochs(EEG_filt1, EEG_filt2, epoch_win)
% Extract epochs

% Required inputs:
%   EEG_fil1 and EEG_filt2: EEG data in eeglab data structure, data is in
% continuous, non-epoched format, containing event timestamps
%   epoch_win: time window relative to event timestamps in seconds, ex: [-1 1];

%-------------------------
[EEG_ep1, indices, com] = pop_epoch(EEG_filt1, {}, epoch_win); % for analyses
%                                             setnm = [curr_resdir filesep 'EEG_ep1_' Lnfilt '_05HP.set'];
%                                             pop_saveset(EEG_ep1,setnm);
if ~isempty(EEG_filt2)
    [EEG_ep2, indices, com] = pop_epoch(EEG_filt2, {}, epoch_win); % for ICA
    
    ep_nr = length(indices);
    
    for i = 1:ep_nr
        EEG_ep2.epoch(i).index = i;
    end
    
else
    EEG_ep2 = [];
end


