function new_Evinxx_parts(curr_resdir,event,newpart_names,newpart_inx)
%NEW_EVINXX_PARTS Adds new fields to Evinxx struct.
% NEW_EVINXX_PARTS(curr_resdir,event,newpart_names,newpart_inx) saves
% indices of new trial sets (NEWPART_INX), intersected with indices of
% epochs aligned to EVENT. Trial/ epoch indiced are stored as
% Evinxx.(EVENT).(NEWPART_NAMES).TE_index/ ...epoch_index.
%
% See also: SAVE_EVINXX
%
% Johanna Petra Szabó, Panna Heged?s, 05.2025
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


part_nr = length(newpart_names);
try
    load(fullfile(curr_resdir,'Evinxx.mat'));
catch
    fprintf('No Evinxx %s.\n',curr_resdir);
    keyboard;
    return;
end

event_TEindex = Evinxx.(event).(event).TE_index;
event_epochindex = Evinxx.(event).(event).epoch_index;

for j = 1:part_nr
    [both_TEindex,einx,~] = intersect(event_TEindex,newpart_inx{j});
    
    Evinxx.(event).(newpart_names{j}).epoch_index = event_epochindex(einx);
    
    Evinxx.(event).(newpart_names{j}).TE_index = both_TEindex;
end



save(fullfile(curr_resdir,'Evinxx.mat'),'Evinxx')