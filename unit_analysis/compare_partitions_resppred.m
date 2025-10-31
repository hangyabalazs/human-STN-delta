function compare_partitions_resppred(EventTypes,partition,parttags,patgr_nm)
%COMPARE_PARTITIONS_RESPPRED Compare partitioned trials in responsive/ predictive units
% COMPARE_PARTITIONS_RESPPRED(EventTypes,partition,parttags)
%   -Compares partitioned sets of trials in responsive/ predictive units 
%       using permutation test with cluster based correction.
%
%Input parameters:
%     EVENTTYPES        1xN cell array of event label for responsive units 
%
%     PARTITION         char. array of partition label (see defineLabelsColors_pd.m)
%
%     PARTTAGS          tag of partition (see defineLabelsColors_pd.m)
%
%
% See also: AVG_PSTH_STAT, DEFINELABELSCOLORS_PD, STD_STAT, NORM_PSTH_MAP1,ULTIMATE_PSTH
%
% Johanna Petra Szabó, 05.2025
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

global cell_dir

resdir = fullfile(cell_dir,'grouped2','Partitions',[partition(2:end)]); 
if ~isdir(resdir); mkdir(resdir); end;

% yL = [-4 4];
yL = [];

for ei = 1:length(EventTypes)
    event = EventTypes{ei};
    avg_psth_stat(event,partition,resdir,[-1.2 3],parttags,patgr_nm,yL)
    
    if strcmp(event,'StopSignal')
        ev2excl = 'StimulusOn';
        avg_psth_stat(event,partition,resdir,[-1.2 3],parttags,patgr_nm,yL,ev2excl)
    end
    
    
end