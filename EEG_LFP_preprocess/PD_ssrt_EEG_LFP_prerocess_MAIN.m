function PD_ssrt_EEG_LFP_prerocess_MAIN(epoch_win,baseline_win, EventTypes,SubEventTypes)
% PD_SSRT_BEHAV_MAIN Main function for EEG and LFP preprocessing 
%
% Input parameters:
%     CONDITIONS      Nx2 cell array of task conditions to compare
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

% See also: PREPROCESS_PD

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

global filesdir


for rT = 1:4
    switch rT
        case 1;
            rectype  = 'EEG';  rectime = 'postop'; condi = 'stimoff'; 
        case 2;
            rectype  = 'EEG';  rectime = 'postop'; condi = 'stimon'; 
        case 3;
            rectype  = 'LFP';  rectime = 'intraop'; condi = 'stimoff'; 
        case 4;
            rectype  = 'EEG';  rectime = 'intraop'; condi = 'stimoff'; 
    end
    
    
          
        sess2analyse = getdata2analyse(filesdir, 'rectype',rectype,...
            'rectime',rectime,'patients', 'allpatients', 'side','bothside', 'condition',condi);
        
        
        preprocess_PD(sess2analyse,epoch_win,EventTypes,SubEventTypes,'epochs',baseline_win);
        
        cuepair_trialtypes_2Evinx(sess2analyse,EventTypes(1));
    end

end