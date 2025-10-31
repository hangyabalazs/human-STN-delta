function spike_phase_coupling_MAIN(EventTypes,SubEventTypes,rectype)
%SPIKE_PHASE_COUPLING
%     Calculates the amount of coupling between spiking activity and LFP/ EEG oscillations (SPC),
%     expressed in mean resultant length (MRL), for each detected unit. SPC is
%     calculated in time windows around a specific behavioral event.
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

%%
% Parameters
rectime = 'intraop'; condi = 'stimoff';

extract_phase_win = [-2 2]; % timestamps of phase values are saved -> smaller time-windows can be retrieved later
plot_win = [-2 0]; % [-2 0] | [0 2] | [-1.5 1.5]; % time period for plotting & analysis

freq_bands = [1 4]; %[1 4];
fr_names1 = {'high_delta'}; % {'beta'} | {'high_delta'};
fr_names2 = {'dom_high_delta'};% {'dom_beta'} | {'dom_high_delta'};

if strcmp(rectype,'LFP')
    group_types = {'signPC', 'all resp', ...
        'StimulusOn resp Active', 'StimulusOn resp Inhib', ...
        'StopSignal resp Active',  'StopSignal resp Inhib',...
        'KeyPress1 resp Active', 'KeyPress1 resp Inhib',...
        'SUA','MUA'};
elseif strcmp(rectype,'EEG')
    group_types = {'signPC', 'signPC_LFP', 'all resp', ...
        'StimulusOn resp Active', 'StimulusOn resp Inhib', ...
        'StopSignal resp Active',  'StopSignal resp Inhib',...
        'KeyPress1 resp Active', 'KeyPress1 resp Inhib',...
        'SUA','MUA'}; % 'RevSkip_slower'

end


% Phase value extraction & Spike-phase coupling measures

sess2analyse = getdata2analyse(filesdir, 'rectype',rectype,...
    'rectime',rectime,'patients', 'allpatients', 'side','bothside', 'condition',condi);


PCdir = fullfile(cell_dir,'PC',[rectype '_phase_coupling']);
if ~isfolder(PCdir); mkdir(PCdir); end;

%%
for s = 1:2
    switch s; case 1; subevs = false; case 2; subevs = true; end;
    
    
    %Phase extraction
    spik_phas_extraction(sess2analyse,extract_phase_win,freq_bands,fr_names1,...
        PCdir,EventTypes,SubEventTypes,'subevs',subevs);
    
    
    
    %Spike-phase coupling measures
    if ~subevs; ds = 'no'; else; ds = 'spike'; end;
    
    PC_cell_level(plot_win,freq_bands, fr_names2,PCdir,EventTypes,SubEventTypes,...
        'downsamp',ds,'subevs',subevs,'rectype',rectype);
    
    
%     % PC around Stop Signal with modified window
%     PC_cell_level([-2 0],freq_bands, fr_names2,PCdir,{'StimulusOn'},{'FailedStopTrial','SuccesfulStopTrial'},...
%         'downsamp',ds,'subevs',subevs,'rectype',rectype);
%     
%     PC_cell_level([0 2],freq_bands, fr_names2,PCdir,{'StimulusOn'},{'FailedStopTrial','SuccesfulStopTrial'},...
%         'downsamp',ds,'subevs',subevs,'rectype',rectype);
%     
    
end

%% Ratios of sign. coupled units
for ei = 1:length(EventTypes)
    event = EventTypes{ei};
    
    % Percentage of sign. coupled units from all detected units (pie charts)
    PC_pie(event,plot_win,fr_names2{1},'no',1,0.05,rectype,false,'all',PCdir,'')
%     PC_pie(event,plot_win,fr_names2{1},'no',1,0.05,rectype,false,'all',PCdir,'SUA')
%    
    PC_pie(event,plot_win,fr_names2{1},'no',1,0.05,rectype,false,'all',PCdir,[event ' resp'])
    
    % Percentage of sign. coupled units within unit groups (stacked bar plots)
    PC_respcells_stacked_compare(event,plot_win,fr_names2{1},'no',1,0.05,rectype,false,'all',PCdir)
    
    
    PC_pie(event,plot_win,fr_names2{1},'no',1,0.05,rectype,false,'all',PCdir,'signPC_dom_beta')
        
        
    if strcmp(rectype,'EEG')
        PC_pie(event,plot_win,fr_names2{1},'no',1,0.05,rectype,false,'all',PCdir,'signPC_LFP')
    end
end


%% Population phase/ mean resultant length distribution in unit groups
for g = 1:length(group_types)
    
    group = group_types{g};
    
    %%
    % Generate plots
    act_rv = PC_groups_f(plot_win,freq_bands,fr_names2,PCdir,EventTypes,SubEventTypes,...
        'rectype',rectype,'group',group,'components',[]);
    
    % Compares subevents (partitions) 
    PC_group_subevs(plot_win,freq_bands,fr_names2,PCdir,EventTypes,SubEventTypes,...
        'group',group,'downsamp','no','partition','#StopPartition','parttags',[1 2]);
    
end

%%
group = 'signPC'; % 'RevSkip_slower' | 'signPC'
PC_groups_f(plot_win,freq_bands,fr_names2,PCdir,EventTypes(1),{'Ord','RevSkip'},...
    'rectype',rectype,'group',group,'components',[1 2 3]);

% Compares subevents (partitions)
PC_group_subevs(plot_win,freq_bands,fr_names2,PCdir,EventTypes(1),{'Ord','RevSkip'},...
    'group',group,'downsamp','no','partition','#CuepairPartition','parttags',[1 4]);
%%
for ei = 1:length(EventTypes)
    event = EventTypes{ei};
    phase_hist_sinewave_plot(event,plot_win,PCdir,15)
end



%% Spike-triggered LFP average of significantly delta-coupled units
spike_triggered_average_PD(sess2analyse,false,fr_names2{1},freq_bands,extract_phase_win,plot_win,EventTypes)



%% RT/SSDp0.5 correaltion with spike-phase coupling
for b = 1:2
    switch b; case 1; param = 'RT'; case 2; param = 'ssd05'; end;
    
    for ei = 1:length(EventTypes)
        event = EventTypes{ei};
        %%
        comp_vals = PC_RT_correlation(param,event, fr_names2{1}, false, rectype,'signPC', true,1,'no',plot_win,'');
%         PC_RT_correlation(param,event, fr_names2{1}, false, rectype,'signPC', true,1,'no',plot_win,'SUA')
%         PC_RT_correlation(param,event, fr_names2{1}, false, rectype,'signPC', true,1,'no',plot_win,'MUA')
    end
    
end

%% Spike-phase coupling and bursting index
PC_bursting(EventTypes,rectype,plot_win)

%%
% Spike-phase coupling and STN localization
PC_STN_loc(EventTypes,rectype, fr_names2{1}, false)


%% PSTH of sign. coupled units
for ei = 1:length(EventTypes)
    event = EventTypes{ei};
    PC_cells_PSTH(event, fr_names2{1}, false, rectype, true,1, 'no', plot_win)
end
end




