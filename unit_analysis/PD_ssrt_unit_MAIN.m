function PD_ssrt_unit_MAIN(EventTypes,EvTy)
% PD_SSRT_BEHAV_MAIN Main function for unit analysis
%
% Input parameters:
%     EVENTTYPES        1xN cell array of event label for responsive units 
%
%     EVTY              1xN cell array of event labels for bursting analyses
%
%
% See also: AUTOCORR_PD, UPDRS_BURSTING_CORR, UNIT_SORTER, UNIT_SUBREGIONS

% Balázs Hangya, Panna Hegedus, Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu




% Single-/ multi units (based on L ratio and ID distribution)
ID_thr = 15;
L_thr = 0.02;
L_ratio_ID_distrib(L_thr,ID_thr)

% Mean firing rate
set_meanFR

% Responsive/ predictive units
unit_sorter(EventTypes,{'StopSignal'})


% Compare partitions of responsive/ predictive units
compare_partitions_resppred(EventTypes,'#StopPartition',[1 2],'')
compare_partitions_resppred(EventTypes,'#StopPartition',[1 2],'psth_stat')
compare_partitions_resppred(EventTypes,'#StopPartition',[1 2],'SUA')
compare_partitions_resppred(EventTypes,'#StopPartition',[1 2],'MUA')

if ismember('StimulusOn',EventTypes)
    patgr_nm = ''; % 'RevSkip_slower' | ''
    compare_partitions_resppred({'StimulusOn'},'#CuepairPartition',[1 4],patgr_nm)
end

% Map of all types of responsive/ predictive units
all_resppred_map(EventTypes,{'StopSignal'})


% Auto-correlation and bursting index
cellids = findcell;
autocorr_PD(cellids,{'delta'},[1 4],3);


% Bursting properties of units
unit_bursting(EvTy,partition,parttags)

% Compare PSTHs of trials separated by median BI
mediansplitBI_PSTHs(EventTypes,'all',[],'')
mediansplitBI_PSTHs(EventTypes,'#StopPartition',[1 2],'')
mediansplitBI_PSTHs({'StimulusOn'},'#StopPartition',[1 2],'SUA');
mediansplitBI_PSTHs({'StimulusOn'},'#StopPartition',[1 2],'MUA');
mediansplitBI_PSTHs({'StimulusOn'},'#CuepairPartition',[1 4],'RevSkip_slower')


% Localization of units within the STN
fish_p = unit_subregions(EventTypes);



end