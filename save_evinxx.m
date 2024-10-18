% Creates structure with indices of good trials following rejection of bad
% trials, according to event types
% + corresponding trial indices (in TrialEvents.mat file)

% EEG_ep: EEG data structure (eeglab format) with event timestamps
% EventTypes: label of events
% SubEventTypes: label of corresponding subevents
% curr_resdir: directory to save Evinxx.mat file
% issave: if true saves Evinxx structure

function Evinxx = save_evinxx(EEG_ep,EventTypes,SubEventTypes,curr_resdir,issave)


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
        if strcmp(event,'StimulusOn') 
        stopei = find(ismember(EventTypes,'StopSignal'));
        for sei = 1:2
            evty = SubEventTypes{stopei,sei};
            
            evinx = [];
            
            [evinx TEevinx] = find_evinx(EEG_ep, evty,3);
            
            if strcmp(evty,'SuccesfulStopTrial') && isempty(evinx)
                [evinx TEevinx] = find_evinx(EEG_ep, 'SuccessfulStopTrial',3);
            end
            
            Evinxx.(event).(evty).epoch_index =evinx;
            Evinxx.(event).(evty).TE_index =TEevinx;
            
        end
    end
    
   [evinx TEevinx] = find_evinx(EEG_ep, event,1);
    Evinxx.(event).(event).epoch_index =evinx;
    Evinxx.(event).(event).TE_index =TEevinx;
end

if issave
    save(fullfile(curr_resdir,['Evinxx.mat']),'Evinxx');
end


