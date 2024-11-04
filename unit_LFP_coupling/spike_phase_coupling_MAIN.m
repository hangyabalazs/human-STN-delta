function spike_phase_coupling_MAIN(EventTypes,SubEventTypes)
%SPIKE_PHASE_COUPLING
%     Calculates the amount of coupling between spiking activity and LFP oscillations (SPC),
%     expressed in mean resultant length (MRL), for each detected unit. SPC is
%     calculated intime windows around a specific behavioral event.
%
%     EVENTTYPES        1xN cell array of event labels, ex: {'StimulusOn','StopSignal'};
%
%     SUBEVENTTYPES     Nx2 cell array of partition ("subevent") labels, each row
%                       corresponds to an event label, columns to partitions,
%                       ex: {'FailedStopTrial','SuccesfulStopTrial';'FailedStopTrial','SuccesfulStopTrial'}
%
% See also: SPIK_PHAS_EXTRACTION, FIND_DOMINANT_FREQ_BANDS, PC_CELL_LEVEL,
%           PC_PIE, PC_RESPCELLS_STACKED_COMPARE,
%           PC_GROUPS, PC_GROUP_SUBEVS, PHASE_HIST_SINEWAVE_PLOT, PC_RT_CORRELATION

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


global filesdir cell_dir

% Parameters


rectype = 'LFP';
rectime = 'intraop'; condi = 'stimoff';

extract_phase_win = [-2 2]; % timestamps of phase values are saved -> smaller time-windows can be retrieved later
plot_win = [-1.5 1.5]; % time period for plotting & analysis

freq_bands = [1 4];
fr_names = {'dom_high_delta'};


group_types = {'signPC', 'all resp', ...
    'StimulusOn resp Active', 'StimulusOn resp Inhib', ...
    'StopSignal resp Active',  'StopSignal resp Inhib'};




% Phase value extraction & Spike-phase coupling measures

sess2analyse = getdata2analyse(filesdir, 'rectype',rectype,...
    'rectime',rectime,'patients', 'allpatients', 'side','bothside', 'condition',condi);


PCdir = fullfile(cell_dir,'PC',[rectype '_phase_coupling']);
if ~isfolder(PCdir); mkdir(PCdir); end;

for s = 1:2
    switch s; case 1; subevs = false; case 2; subevs = true; end;
    
    
    %Phase extraction
    spik_phas_extraction(sess2analyse,extract_phase_win,freq_bands,fr_names,...
        PCdir,EventTypes,SubEventTypes,'subevs',subevs);
    
    
    
    %Spike-phase coupling measures
    if ~subevs; ds = 'no'; else; ds = 'spike'; end;
    
    PC_cell_level(plot_win,freq_bands, fr_names,PCdir,EventTypes,SubEventTypes,...
        'downsamp',ds,'subevs',subevs,'rectype',rectype);
    
end

% Ratios of sign. coupled units
for ei = 1:length(EventTypes)
    event = EventTypes{ei};
    
    % Percentage of sign. coupled units from all detected units (pie charts)
    PC_pie(event,plot_win,fr_names{1},'no',1,0.05,rectype,false,'all',PCdir)
    
    % Percentage of sign. coupled units within unit groups (stacked bar plots)
    PC_respcells_stacked_compare(event,plot_win,fr_names{1},'no',1,0.05,rectype,false,'all',PCdir)
    
end


%Population phase/ mean resultant length distribution in unit groups
for g = 1:length(group_types)
    
    group = group_types{g};
    
    
    % Generate plots
    PC_groups_f(plot_win,freq_bands,fr_names,PCdir,EventTypes,SubEventTypes,...
        'rectype',rectype,'group',group);
    
    % Compares subevents (partitions) based on Stop trial outcome
    PC_group_subevs(plot_win,freq_bands,fr_names,PCdir,EventTypes,SubEventTypes,...
        'group',group,'downsamp','no');
    
end

for ei = 1:length(EventTypes)
    event = EventTypes{ei};
    phase_hist_sinewave_plot(event,plot_win,PCdir,15)
end



% Spike-triggered LFP average of significantly delta-coupled units
spike_triggered_average_PD(sess2analyse,false,freq_bands,extract_phase_win,plot_win,EventTypes)



% RT/SSDp0.5 correaltion with spike-phase coupling
for b = 1:2
    switch b; case 1; param = 'RT'; case 2; param = 'ssd05'; end;
    
    for ei = 1:length(EventTypes)
        event = EventTypes{ei};
        PC_RT_correlation(param,event, fr_names{1}, false, rectype,'signPC', true,1,'no',plot_win)
    end
    
end

% Spike-phase coupling and bursting index
PC_bursting(EventTypes)


end




