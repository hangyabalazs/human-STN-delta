function [evinx TEevinx] = find_evinx(EEG_ep, evty,evlevel)
% FIND_EVINX Gets epoch index corresponding to EVTY event and EVLEVEL level
% of subevents in EEG_EP data (eeglab struct). 
%
%
% Johanna Petra Szabó, 10.2025
% Clinic for Neurosurgery and Neurointervention
% Semmelweis University, Budapest, Hungary
% szabo.johanna.petra@semmelweis.hu

nargoutchk(1,2)

evinx = []; TEevinx = [];
for ive = 1:size(EEG_ep.epoch,2)
    if evlevel==1
        epev = [EEG_ep.epoch(ive).eventtype];
    elseif evlevel==2
        
        epev = [EEG_ep.epoch(ive).eventsubtype];
    elseif evlevel==3
        try
            epev = [EEG_ep.epoch(ive).eventstoptype];
        catch
            keyboard;
        end
    end
    if isempty(epev); continue; end;
    
    if iscell(EEG_ep.epoch(ive).eventlatency);
        central_evix = find([EEG_ep.epoch(ive).eventlatency{:}] == 0);
        TEi = EEG_ep.epoch(ive).eventTEindex{1,central_evix};
        okeix = any(strcmp(evty,epev(central_evix)));
        
    else
        if EEG_ep.epoch(ive).eventlatency==0;
            okeix = strcmp(evty,epev);
            TEi = EEG_ep.epoch(ive).eventTEindex;
        end;
    end
    if okeix
        evinx = [evinx ive];
        TEevinx = [TEevinx TEi];
        
    end
end