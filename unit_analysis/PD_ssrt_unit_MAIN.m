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

% Map of all types of responsive/ predictive units
all_resppred_map(EventTypes,{'StopSignal'})


% Auto-correlation and bursting index
cellids = findcell;
autocorr_PD(cellids,{'delta'},[1 4],3);


% Bursting properties of units
unit_bursting(EvTy)


% Localization of units within the STN
fish_p = unit_subregions(EventTypes)



end