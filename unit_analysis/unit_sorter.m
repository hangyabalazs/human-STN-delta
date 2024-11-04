function unit_sorter(EventTypes,EvTy)
%UNIT_SORTER Responsive/ predictive units
% UNIT_SORTER(EventTypes,EvTy)
%   -Tests activation/ inhibition of units around events
%   -Statistically significant activation/ inhibition following an event
%   (responsive units)
%   or significant difference between partitioned sets of trials before an
%   event (predictive units) are marked. 
%   -Generates average PSTHs and PSTH maps of responsive/ predictive units.
%   -Compares partitioned sets of trials in responsive/ predictive units.
%
%Input parameters:
%     EVENTTYPES        1xN cell array of event label for responsive units 
%
%     EVTY              1xN cell array of event labels for predictive units
%
% See also: RESPONSESORTER_PD, EVSI_AVG_PSTH, RESPSORT_PARTITIONS_PD,
% PARTITIONS_AVG_PSTH, AVG_PSTH_STAT


% Balázs Hangya, Panna Heged?s, Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu



global cell_dir stat_dir group_dir




test_window = [0 1];
twinds = [-0.5 0; -1 0; -1.5 0; -2 0];
pdcells = findcell;


% Find responsive cells
for iE = 1:length(EventTypes)
    
    event = EventTypes{iE};

    responsesorter_PD(pdcells, stat_dir, event, test_window, true);
    
    
end


% Make PSTH maps, PSTH average plots of grouped responsive units

for eii = 1:length(EventTypes)
    event = EventTypes{eii};
    
    propname_resp = [event 'psth_stat_' num2str(test_window)];
    
    EvsI_avg_psth(pdcells,propname_resp,event,'none',true)
    
end



% Find predictive cells



for i = 1:length(EvTy)
    alignevent = EvTy{i};
    resdir = fullfile(cell_dir,'Predict_PD_multi');
    
    respsort_partitions_PD(pdcells,alignevent,'#StopPartition', twinds, resdir, 1);
    
end


% Make PSTH maps, PSTH average plots of grouped predictive units


for i =1:length(EvTy)
    
    alignevent = EvTy{i};
    fprintf('%s...',alignevent)
    propname_resp = [alignevent 'psth_stat_' num2str(test_window)];
        
    partitions_avg_psth(pdcells,twinds,alignevent,'#StopPartition',propname_resp,{'none'},1);
    
end


% Compare partitions of responsive/ predictive units
partition = '#StopPartition';
resdir = fullfile(cell_dir,'grouped2','Partitions',[partition(2:end) '2']);
if ~isdir(resdir); mkdir(resdir); end;


for ei = 1:length(EventTypes)
    avg_psth_stat(EventTypes{ei},'#StopPartition',resdir,[-1.5 3])
end


end