function PD_SSRT_MAIN(rootdir, filesdir)
% PD_SSRT_BEHAV_MAIN Main function for analysis of Stop Signal
% Reaction Time task of Parkinsonian patients implanted with deep brain
% stimulator. 
%
% Input parameters:
%       ROOTDIR       path to directory to save results
%
%       FILESDIR      path to directory containing files to analyse 

% Johanna Petra Szabó, Panna Hegedus, Balázs Hangya, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu



% Directories
figdir_pd = fullfile(rootdir,'Figures_epoched'); if ~isfolder(figdir_pd); mkdir(figdir_pd); end;
cell_dir = fullfile(figdir_pd,'Intraop_SP');  if ~isfolder(cell_dir); mkdir(cell_dir); end;
group_dir = fullfile(cell_dir,'grouped2\psth_stat1');  if ~isfolder(group_dir); mkdir(group_dir); end;
stat_dir = fullfile(cell_dir,'Response_PD','PD_responsesorter1');  if ~isfolder(stat_dir); mkdir(stat_dir); end;

global rootdir filesdir figdir_pd cell_dir group_dir stat_dir

% Parameters
epoch_win = [-2 2];
baseline_win = [-1 -.5];

EventTypes_all= {'StimulusOn','StopSignal','KeyPress1','Feedback'};
SubEventTypes_all = {'CueStim','StopStim';'FailedStopTrial','SuccesfulStopTrial';'CueResponse','StopResponse';'Correct','Error'};

EventTypes = {'StimulusOn','StopSignal'};
SubEventTypes = {'FailedStopTrial','SuccesfulStopTrial';'FailedStopTrial','SuccesfulStopTrial';};


task_conditions = {'preop','stimoff';'intraop','stimoff';'postop','stimoff';'postop','stimon'};
all_groups =  {{'all'}, {'tremor-dominant','akinetic-rigid','mixed'} , {'RTdecrease','RTincrease'}};


% Behavioural analysis 
PD_ssrt_behav_MAIN(task_conditions, all_groups)

% Preprocess EEG and LFP data
PD_ssrt_EEG_LFP_prerocess_MAIN(task_conditions(2:end,:), epoch_win,baseline_win, EventTypes_all,SubEventTypes_all)

% Time-frequency analysis of EEG and LFP data
PD_ssrt_EEG_LFP_time_freq_MAIN(task_conditions(2:end,:), epoch_win,baseline_win, EventTypes,SubEventTypes)

% Intraop EEG - LFP wavelet coherence analysis
PD_EEG_LFP_wav_coherence_MAIN(EventTypes, SubEventTypes)

% Unit analysis
PD_ssrt_unit_MAIN(EventTypes_all,EventTypes)

% Coupling of unit spiking activity with LFP delta oscillation
spike_phase_coupling_MAIN(EventTypes,SubEventTypes)

% Comparison of patient groups (patients with reaction time increase and -decrease
patient_groups_compare(EventTypes)

