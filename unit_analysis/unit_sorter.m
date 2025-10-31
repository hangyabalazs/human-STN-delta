function unit_sorter(EventTypes,EvTy)
%UNIT_SORTER Responsive/ predictive units
% UNIT_SORTER(EventTypes,EvTy)
%   -Tests activation/ inhibition of units around events
%   -Statistically significant activation/ inhibition following an event
%   (responsive units)
%   or significant difference between partitioned sets of trials before an
%   event (predictive units) are marked. 
%   -Generates average PSTHs and PSTH maps of responsive/ predictive units.
%
%Input parameters:
%     EVENTTYPES        1xN cell array of event label for responsive units 
%
%     EVTY              1xN cell array of event labels for predictive units
%
% See also: RESPONSESORTER_PD, EVSI_AVG_PSTH, RESPSORT_PARTITIONS_PD, 
% PARTITIONS_AVG_PSTH, AVG_PSTH_STATNORM_PSTH_MAP1,ULTIMATE_PSTH


% Balázs Hangya, Panna Heged?s, Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu



global cell_dir stat_dir group_dir




test_window = [0 1];
twinds = [-0.5 0; -1 0; -1.5 0; -2 0];
pdcells = findcell;
isfig = true;
% Find responsive cells
for iE = 1:length(EventTypes)
    
    event = EventTypes{iE};
    
    
    sort_resdir = fullfile(stat_dir,event); if ~isdir(sort_resdir); mkdir(sort_resdir); end;
    responsesorter_PD(pdcells, sort_resdir, event, test_window, true,'',isfig);
    
    if strcmp(event,'StopSignal')
        ev2excl = 'StimulusOn';
        sort_resdir = fullfile(stat_dir,[event '_ExclBef_' ev2excl]); if ~isdir(sort_resdir); mkdir(sort_resdir); end;
        responsesorter_PD(pdcells, sort_resdir, event, test_window, true,ev2excl,isfig);
    end
    
    
end


% Make PSTH maps, PSTH average plots of grouped responsive units

for eii = 1:length(EventTypes)
    event = EventTypes{eii};
    
    propname_resp = [event 'psth_stat_' num2str(test_window)];
    propname_boxstat = [event 'boxstat_' num2str(test_window)];
    
%     EvsI_avg_psth(pdcells,{propname_resp, propname_boxstat},event,'none',true)
    EvsI_avg_psth(pdcells,propname_resp,event,'none',true,[],true,'_psth_stat')
    
    if strcmp(event,'StopSignal')
        ev2excl = 'StimulusOn';
        EvsI_avg_psth(pdcells,{propname_resp, propname_boxstat},event,'none',true,ev2excl)
    end
    
    
end



% Find predictive cells



for i = 1:length(EvTy)
    alignevent = EvTy{i};
    resdir = fullfile(cell_dir,'Predict_PD_multi');
    
    respsort_partitions_PD(pdcells,alignevent,'#StopPartition', twinds, resdir, 1,[]);
    
    if strcmp(alignevent,'StopSignal')
        ev2excl = 'StimulusOn';
        
        respsort_partitions_PD(pdcells,alignevent,'#StopPartition', twinds, resdir, 1,ev2excl);
    end
    
end


% Make PSTH maps, PSTH average plots of grouped predictive units


for i =1:length(EvTy)
    
    alignevent = EvTy{i};
    fprintf('%s...',alignevent)
    propname_resp = [alignevent 'psth_stat_' num2str(test_window)];
    
    partitions_avg_psth(pdcells,twinds,alignevent,'#StopPartition',propname_resp,{'none'},1);
    
    if strcmp(alignevent,'StopSignal')
        ev2excl = 'StimulusOn';
        partitions_avg_psth(pdcells,twinds,alignevent,'#StopPartition',propname_resp,{'none'},1,ev2excl);
    end
end



end