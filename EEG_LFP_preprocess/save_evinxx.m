function Evinxx = save_evinxx(EEG_ep,EventTypes,SubEventTypes,curr_resdir,issave)
%SAVE_EVINXX    Saves good trial indeces
%   Evinxx = SAVE_EVINXX(EEG_ep,EventTypes,SubEventTypes,curr_resdir,issave)
%       creates struct with fields corresponding to events (EVENTTYPES) and 
%       subevents (SUBEVENTTYPES). The subfield 'epoch_index' contains indeces 
%       of event epochs in EEG_EP.DATA (3D matrix: channels x time x epochs).
%       The subfield 'TE_index' contains indeces of the corresponding
%       events in TrialEvents.mat file (loaded from the result directory of
%       the patient: CURRSESS). 
%       Only good epochs remaining following rejection of bad trials are
%       included. 
%       The resulting struct (EVINXX) is saved as Evinxx.mat file in the
%       result directory of the patient (CURR_RESDIR), if ISSAVE is true.
% 
% See also: PREPROCESS_PD

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


for ei = 1:length(EventTypes)
    event = EventTypes{ei};
    
    for sei = 1:2
        evty = SubEventTypes{ei,sei};
        
        evinx = [];
        
        [evinx TEevinx] = find_evinx(EEG_ep, evty,2);
        
        Evinxx.(event).(evty).epoch_index =evinx;
        Evinxx.(event).(evty).TE_index =TEevinx;
                
    end
    
%     if strcmp(event,'StimulusOn') || strcmp(event,'KeyPress1')
%         if strcmp(event,'StimulusOn') 
%         stopei = find(ismember(EventTypes,'StopSignal'));
%         for sei = 1:2
%             evty = SubEventTypes{stopei,sei};
%             
%             evinx = [];
%             
%             [evinx TEevinx] = find_evinx(EEG_ep, evty,3);
%             
%             if strcmp(evty,'SuccesfulStopTrial') && isempty(evinx)
%                 [evinx TEevinx] = find_evinx(EEG_ep, 'SuccessfulStopTrial',3);
%             end
%             
%             Evinxx.(event).(evty).epoch_index =evinx;
%             Evinxx.(event).(evty).TE_index =TEevinx;
%             
%         end
%     end
    
   [evinx TEevinx] = find_evinx(EEG_ep, event,1);
    Evinxx.(event).(event).epoch_index =evinx;
    Evinxx.(event).(event).TE_index =TEevinx;
end

if issave
    save(fullfile(curr_resdir,['Evinxx.mat']),'Evinxx');
end


