function [evinx TEevinx] = find_evinx(EEG_ep, evty,evlevel)

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
    if any(strcmp(evty,epev([EEG_ep.epoch(ive).eventlatency{:}] == 0)))
        evinx = [evinx ive];
        
        try
            zeroev = find([EEG_ep.epoch(ive).eventlatency{:}] == 0);
            TEi = EEG_ep.epoch(ive).eventTEindex{1,zeroev};
            TEevinx = [TEevinx TEi];
        catch
            if ive==1
                fprintf('No TEindex field %s\n',EEG_ep.filepath);
            end
            continue
            
        end
    end
end