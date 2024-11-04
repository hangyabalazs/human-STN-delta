function create_EventSpikes_mat(pdcells)
%CREATE_EVENTSPIKES_MAT     Saves event-related spikes and epoch rates
%   CREATE_EVENTSPIKES_MAT(pdcells) 
%       -define events and event epochs for SSRT task
%       -extract and align spikes in each trial relative to trial events
%       -calculate spike rates in fixed or variable epochs.
%       -save event spikes and epoch rates in a struct as
%       'EVENTSPIKES*.MAT' in the session's result directory.
%
% See also: DEFINEEVENTSEPOCHS_PDTASK, EXTRACTEVENTSPIKES, EXTRACTEPOCHRATES

% Balázs Hangya, Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary


for ic = 1:length(pdcells)
    cellid = pdcells{ic};
   TE = loadcb(cellid,'TrialEvents');
    [events,epochs] = defineEventsEpochs_pdtask; 
    
    [event_stimes,event_windows] = extractEventSpikes(cellid,events,TE);
    
    [epoch_rates,~,~] = extractEpochRates(event_stimes,event_windows,events,epochs);
    
    EvSp.cellid = cellid;
    EvSp.epoch_rates = epoch_rates;
    EvSp.epochs = epochs;
    EvSp.event_stimes = event_stimes;
    EvSp.event_windows = event_windows;
    EvSp.events = events;
    
    cellbase_datapath = getpref('cellbase','datapath');
    [patname,session,tetrode,unit] = cellid2tags(cellid);
    
    cell_fname = fullfile(cellbase_datapath,patname,session);
    
    save(fullfile(cell_fname,['EVENTSPIKES' num2str(tetrode) '_' num2str(unit)]),'-struct','EvSp');
    
end