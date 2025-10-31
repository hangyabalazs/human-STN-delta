function [EEGevs] = behav_events(EEG,EventTypes,SubEventTypes,currsess,tag)
%BEHAV_EVENTS   Saves event timestamps to EEG (eeglab) data structure
% [EEGevs] = BEHAV_EVENTS(EEG,EventTypes,SubEventTypes,currsess,tag) 
%       loads TrialEvents.mat (contains timestamps already synchronized 
%       between EEG system and task computer)  from the patient's result
%       directory (CURRSESS), stores event (EVENTTYPES) and subevent
%       (SUBEVENTTYPES) timestamps in 'event' field of EEGevs data structure. 
%
% Input parameters: 
%   EEG             eeglab data structure, contains continuous EEG data
% 
%   EVENTTYPES      label of events (same as the field names in TrialEvents.mat)
% 
%   SUBEVENTTYPES   label of corresponding subevents (as in TrialEvents.mat)
% 
%   CURRSESS        patient result directory, contains TrialEvents.mat file
%   
%   TAG             task condition ('stimoff' | 'stimon')
%   
% Output parameters:
%   EEGEVS          eeglab data structure, with 'event' filed populated
%                   with event information

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


load([currsess filesep 'TrialEvents_' tag '.mat']);


Events = struct;
EEGevs = EEG; EEGevs.origevent = EEG.event; EEGevs.event = [];

for i = 1:length(EventTypes)
    Events.(EventTypes{i}) = (TE.TrialStart + TE.(EventTypes{i}))*EEG.srate;
    
    
    for k = 1:length(Events.(EventTypes{i})) % loop over trials
        
        if isnan(Events.(EventTypes{i})(k))
            continue
        end
        try
            subevinx  = find(~isnan([TE.(SubEventTypes{i,1})(k) TE.(SubEventTypes{i,2})(k)]));
            if isempty(subevinx) && strcmp(EventTypes{i},'Feedback')
                fprintf('%s, %d',EventTypes{i},k);
                if TE.FailedStopTrial(k)==1
                    TE.Error(k)=1;
                elseif TE.SuccessfulStopTrial(k)==1
                    TE.Correct(k)=1;
                end
                subevinx  = find(~isnan([TE.(SubEventTypes{i,1})(k) TE.(SubEventTypes{i,2})(k)]));
            end
            
            if strcmp(EventTypes{i},'StimulusOn')
                if ~isnan(TE.StopSignal(k))
                    if TE.SuccesfulStopTrial(k)==1
                        stoptype = 'SuccesfulStopTrial';
                    elseif TE.FailedStopTrial(k)==1
                        stoptype = 'FailedStopTrial';
                    else
                        stoptype = 'NoStopTrial';
                    end
                else
                    stoptype = 'NoStopTrial';
                end
            else
                stoptype = '';
            end
            
            
            if isempty(EEGevs.event)
                EEGevs.event(1).latency = Events.(EventTypes{i})(k);
                EEGevs.event(1).type = EventTypes{i};
                EEGevs.event(1).subtype = SubEventTypes{i,subevinx};
                EEGevs.event(1).index = 1;
                EEGevs.event(1).stoptype = stoptype;
            else
                EEGevs.event(end+1).latency =  Events.(EventTypes{i})(k);
                EEGevs.event(end).type = EventTypes{i};
                EEGevs.event(end).subtype = SubEventTypes{i,subevinx};
                EEGevs.event(end).index = length(EEGevs.event);
                EEGevs.event(end).stoptype = stoptype;
            end
            
        catch
            disp(k);
        end
        
        EEGevs.event(end).TEindex = k;
        
    end
end