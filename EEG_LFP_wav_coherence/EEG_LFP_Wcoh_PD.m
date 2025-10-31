function EEG_LFP_Wcoh_PD(EventTypes,side)
% EEG_LFP_WCOH_PD  Mean-squared wavelet coherence (MSWC) map between intraop EEG and LFP
%   EEG_LFP_WCOH_PD(EventTypes,SubEventTypes) draws an MSWC map
%   patient-by-patient and averaged across patients around EVENTTYPES.
%   Compares epochs partitioned based on SUBEVENTTYPES.
%
%   Input parameters:
%     EVENTTYPES        1xN cell array of event labels, ex: {'StimulusOn','StopSignal'};
%
%     SUBEVENTTYPES     Nx2 cell array of partition ("subevent") labels, each row
%                       corresponds to an event label, columns to partitions
%                       {'FailedStopTrial','SuccesfulStopTrial';'FailedStopTrial','SuccesfulStopTrial'};
%
% See also: WCOHERENCE, STATCONDFIELDTRIP

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


global filesdir

% Parameters
plot_win = [-1 1];
baseline_win = [-1 -.5];
all_freqs = [1 80];


sess2analyse = getdata2analyse(filesdir, 'rectype','EEG',...
    'rectime','intraop','patients', 'allpatients', 'side',side, 'condition','stimoff');



%% WCoherence patient-by-patient - average across event epochs (plot + save struct)
for ei = 1:length(EventTypes)
    
    event = EventTypes{ei}; fprintf('Event: %s,',event);
    
    
    subevent = '';
    downsamp = false;
    
    wcoh_onebyone(sess2analyse,event,subevent,0,all_freqs,downsamp,NaN,plot_win,[],1);
%     wcoh_onebyone(sess2analyse,event,subevent,1,all_freqs,downsamp,0.05,plot_win,baseline_win,1);
    
end



%% WCoherence - average across patients (plot + save struct)
for ei = 1:length(EventTypes)
    
    event = EventTypes{ei}; fprintf('Event: %s,',event);
    
    wcoh_avg(event,'freqs',all_freqs,'alpha',0.05,'plot_win',plot_win,'baseline_win',baseline_win,'side',side);
    
end

% Get baseline coherence value (not normalized)
wcoh_avg(event,'freqs',all_freqs,'alpha',NaN,'plot_win',plot_win,'baseline_win',[],'side',side);

%% Time series avg. across delta band/ subbands
for ei = 1:length(EventTypes)
    
    event = EventTypes{ei}; fprintf('Event: %s,',event);
    wcoh_avg(event,'freqs',[1 4],'alpha',NaN,'plot_win',plot_win,'baseline_win',baseline_win,'side',side,'avgband',true);
    
end

end






%--------------------------------------------------------------------------






