function [EEG_ep1,EEG_ep2] =  prep_epochs(EEG_filt1, EEG_filt2, epoch_win)
%PREP_EPOCHS        Prepares data epochs
%   [EEG_ep1,EEG_ep2] =  prep_epochs(EEG_filt1, EEG_filt2, epoch_win)
%   extract data epochs from continuous EEG data (EEG_FILT1, EEG_FILT2)
%   stored as eeglab data structure, using event information stored in
%   the 'event' field of these structs. Epochs include the time windows
%   defined by EPOCH_WIN (sec) around each event. Epoched data is stored 
%   in 'data' field of structs as 3D matrix: channels x time x epochs.
% 
% 
%   Input parameters:
%       EEG_FILT1, EEG_FILT2        structs, eeglab data structure; 
%                                   'data' field contains 2D matrix (channels x time) of
%                                   continuous, non-epoched recording
%                                   'event' field contains event timestamps
%
%       EPOCH_WIN:                  time window relative to event timestamps in seconds, ex: [-1 1];
%
%   Output parameters:
%       EEG_EP1, EEG_EP2            same format as EEG_FILT1, EEG_FILT2 structs with epoched data 
%
% See also: PREPROCESS_PD

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


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


