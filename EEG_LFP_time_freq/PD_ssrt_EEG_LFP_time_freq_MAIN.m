function PD_ssrt_EEG_LFP_time_freq_MAIN(conditions, epoch_win,baseline_win, EventTypes,SubEventTypes)
% PD_SSRT_BEHAV_MAIN Main function for EEG and LFP time-frequency analysis
%
% Input parameters:
%       CONDITIONS      Nx2 cell array of task conditions to compare
%               -first column corresponds to recording time, second column to
%               DBS stimulation condition
%               -ex.:{'preop','stimoff';'intraop','stimoff';'postop','stimoff';'postop','stimon'};
%     EPOCH_WIN         1x2 vector, time window relative to event timestamp in sec, for data epoching (ex: [-2 2])
%
%     BASELINE_WIN       1x2 vector, time window relative to event timestamp in sec, for baseline correction
%
%     EVENTTYPES        1xN cell array of event labels, ex: {'StimulusOn','StopSignal','KeyPress1','Feedback'};
%
%     SUBEVENTTYPES     Nx2 cell array of partition ("subevent") labels, each row
%                       corresponds to an event label, columns to partitions
%                       {'CueStim','StopStim';'FailedStopTrial','SuccesfulStopTrial';'CueResponse','StopResponse';'Correct','Error';};

% See also: TIME_FREQ_PATIENTS, ERSP_PLOT_STAT, PD_EEG_STATS,
%           FIND_DOMINANT_FREQ_BANDS, TFPOWER_MAP_RT

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


global filesdir

% Parameters
condnr = size(conditions,1);
alphas = 0.05;
choi = {};
stat_time = [-1 1];

% Parameters for ERPS plots
subregion = 'all'; % 'all' | 'Motor' | 'Limbic' | 'Associative'
chanmean = 1;

topo_freq_nms = {'delta'};
topo_freqs = [1 4];
topobin = 500;

for rT = 1:3
    switch rT
        case 1;
            rectype  = 'EEG';  rectime = 'postop';
            csd = true;    bipol = false; condi = 'bothcond'; cnr = 2;
        case 2;
            rectype  = 'EEG';
            csd = false; bipol = true;
            rectime = 'intraop'; condi = 'stimoff'; cnr = 1;
        case 3;
            rectype = 'LFP';
            csd = false; bipol = false;
            rectime = 'intraop'; condi = 'stimoff'; cnr = 1;
    end
    
    
    sess2analyse = getdata2analyse(filesdir, 'rectype',rectype,...
        'rectime',rectime,'patients', 'allpatients', 'side','bothside', 'condition',condi);
    
    
    % Time-frequency decomposition
    
    time_freq_patients(sess2analyse,EventTypes,SubEventTypes,epoch_win,baseline_win,...
        'stat_time',stat_time,'csd',csd,'bipol',bipol,'alpha',alphas)
    
    % ERSP plots
    for si = 1:2
        switch si; case 1; side = 'left'; case 2; side = 'right'; end;
        
        
        
        for ci = 1:cnr
            switch ci; case 1; tag = 'stimoff'; case 2; tag = 'stimon'; end
                  
            sess2analyse_s = sess2analyse(ismember({sess2analyse.side},side));
            sess2analyse_sc = sess2analyse_s(ismember({sess2analyse_s.tag},tag));
            
            ERSP_plot_stat(sess2analyse_sc,EventTypes,SubEventTypes, epoch_win,baseline_win,...
                'subregion',subregion,'csd',csd,'bipol',bipol,'alpha',alphas,'stat_time',stat_time,...
                'side',side,'condition',tag,'chanmean',chanmean,...
                'topo_freq_nms',topo_freq_nms,'topo_freqs',topo_freqs,'topobin',topobin);
            
        end
        
        
    end
    
    
    
    % STAT
    if rT==1
        groups2comp_labels = {'conditions','partitions'};
    else
        groups2comp_labels = {'partitions'};
    end
    
    for si = 1:2
        switch si; case 1; side = 'left'; case 2; side = 'right'; end;
        sess2analyse_s = sess2analyse(ismember({sess2analyse.side},side));
        
        
        PD_eeg_stats(sess2analyse_s,EventTypes,SubEventTypes,groups2comp_labels,...
            'side',side,'csd',csd,'bipol',bipol,'fr_names',topo_freq_nms, 'freq_bands',topo_freqs);
    end
    
    
    
    % Find dominant frequency band within delta range
    find_dominant_freq_bands(sess2analyse,topo_freq_nms,topo_freqs,chanmean,subregion);
    
    
    
end

% Correlation map between reaction time and wavelet power coefficients
TFpower_map_RT

end