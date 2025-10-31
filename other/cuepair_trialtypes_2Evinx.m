function cuepair_trialtypes_2Evinx(sess2analyse,EventTypes)
%CUEPAIR_TRIALTYPES_2EVINX saves trial/epoch indeces of trials with different level of conflict.
% CUEPAIR_TRIALTYPES_2EVINX sorts trials into following groups based on cue pair:
%
% - ordered, i.e. numbers are in ascendent order: 12, 23 (level of
% conflict: lowest), name of trial set: 'Ord'
% - ordered but one number is skipped: 13 (level of conflict: middle),
% name of trial set: 'OrdSkip'
% - reversed order, i.e. numbers are in descendent order: 32, 21 (level of
% conflict: middle), name of trial set: 'Rev'
% - reversed order and one number is skipped: 31 (level of conflict:
% highest), name of trial set: 'RevSkip'
%
% Indeces of sorted trials are saved in Evinxx.mat file (see save_evinxx.m).
%
% Input params:
%     SESS2ANALYSE      struct containing all necessary information (name of patient, side
%                       of experiment, tag of condition, session folder path) of
%                       session data that need to be analysed (see getdata2analyse)
%
%     EVENTTYPES        1xN cell array of event labels, ex: {'StimulusOn','StopSignal','KeyPress1','Feedback'};
%
% See also: CUEPAIR_TRIALTYPES, NEW_EVINXX_PARTS, SAVE_EVINXX
%
% Johanna Petra Szabó, Panna Heged?s, 05.2025
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna.petra@semmelweis.hu

for i = 1:length(sess2analyse)
    
    sessdir = sess2analyse(i).sessfolder;
    curr_resdir = sess2analyse(i).folder;
    condition = sess2analyse(i).tag;
    
    TEnm = ['TrialEvents_' condition '.mat'];
    
    % Get trialtypes based on presented cue pairs
    try
        [ord, all_ord, ord_skip, rev, rev_skip] = cuepair_trialtypes(sessdir,TEnm);
    catch
        fprintf('No TE % s %s %s\n',sess2analyse(i).patient,...
            sess2analyse(i).side,sess2analyse(i).tag);
        continue;
    end
    
    
    
    for ei = 1:length(EventTypes)
        event = EventTypes{ei};
        
        new_Evinxx_parts(curr_resdir,event,...
            {'Ord','AllOrd','OrdSkip','Rev','RevSkip'},{ord, all_ord,ord_skip,rev, rev_skip});
        
    end
    
    
    TE = load(fullfile(sessdir,TEnm));
    
    if strcmp(fieldnames(TE),{'TE'})
        TE = TE.TE;
    end
    
    [TE.Ord, TE.AllOrd, TE.OrdSkip, TE.Rev, TE.RevSkip] = deal( nan(1, length( TE.StimulusOn )) );
    
    TE.Ord(ord) = 1;
    TE.AllOrd(all_ord) = 1;
    TE.OrdSkip(ord_skip) = 1;
    TE.Rev(rev) = 1;
    TE.RevSkip(rev_skip) = 1;
    
    save(fullfile(sessdir,TEnm),'TE');
    
    
    
    if strcmp(sess2analyse(1).rectime,'intraop') && strcmp(sess2analyse(1).rectype,'LFP')
        
        sessdir = sess2analyse(i).sessfolder;
        TE = load(fullfile(sessdir,'TrialEvents.mat'));
        
        trnr = length(TE.StimulusOn);
        TE.CuepairPartition = nan(1,trnr);
        TE.CuepairPartition(ord) = 1; % Ord
        TE.CuepairPartition(ord_skip) = 2; % OrdSkip
        TE.CuepairPartition(rev) = 3; % Rev
        TE.CuepairPartition(rev_skip) = 4; % RevSkip
%         TE.CuepairPartition(all_ord) = 11; % All Ord
        
        save(fullfile(sessdir,'TrialEvents.mat'),'-struct','TE');
        
    end
    
    
    
    
end


