function PD_EEG_LFP_wav_coherence_MAIN(EventTypes, SubEventTypes)
% PD_EEG_LFP_WAV_COHERENCE_MAIN Main function for intraop EEG - LFP wavelet coherence analysis
% Input parameters
%     EVENTTYPES        1xN cell array of event labels, ex: {'StimulusOn','StopSignal'};
%
%     SUBEVENTTYPES     Nx2 cell array of partition ("subevent") labels, each row
%                       corresponds to an event label, columns to partitions,
%                       ex: {'FailedStopTrial','SuccesfulStopTrial';'FailedStopTrial','SuccesfulStopTrial'}
%
% See also: EEG_LFP_ETA_PD, EEG_LFP_WCOH_PD, WCOH_MAP_RT

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

% Intraop EEG and LFP event-triggered averages (ETA) aligned
% EEG_LFP_ETA_PD(EventTypes)

% Mean-squared wavelet coherence (MSWC) map between intraop EEG and LFP

side = 'left';
EEG_LFP_Wcoh_PD(EventTypes,side)

% Compare partitions
EEG_LFP_Wcoh_compare_partitions(EventTypes,SubEventTypes,side,'#StopPartition','all','common')

if ismember('StimulusOn',EventTypes)
    patgr_nm = 'all'; % 'RevSkip_slower' | 'all'
    EEG_LFP_Wcoh_compare_partitions({'StimulusOn'},{'Ord','RevSkip'},side,'#CuepairPartition',patgr_nm,'common');
end

% Correlation map between reaction time and wavelet coherence map
wcoh_map_RT
end