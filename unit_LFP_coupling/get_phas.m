function [PC_results_dir, mrl, phas, pdcellids,frnm,varargout] = get_phas(event, fr_name, chanmean, rectype,group, dominantfreq,PC_win, ds, plot_win)
%GET_PHAS   Phase coupling parameters
% GET_PHAS(...) loads results of spike phase coupling analysis.
%     
% Input:
%
% EVENTTYPES        char. array of event label, ex: 'StimulusOn','StopSignal','KeyPress1','Feedback'
% 
% FR_NAME         char. array, label of frequency bands of interest
% 
% CHANMEAN     if true, channel averaged LFP data is used
% if false, channel-by-channel LFP data is considered
%     
% RECTYPE       recording type to use for analyses ('EEG' | 'LFP')
% 
% GROUP        char. array of unit groups;
% 'all' | 'signPC' | 'signPC_LFP' ('signPC_LFP' relevant in case of EEG PD analysis, selects unit sign. coupled to LFP)
% [EVENT ' resp'] | [EVENT ' pred'], (EVENT is a char. array of an event label)
% 'UA' (loads single- and multi-units separately)
% 
% DOMINANTFREQ  if true, dominant frequency within predefined frequency
% band limits (FREQS) are used for PC calculation, see FIND_DOMINANT_FREQ_BANDS
% true | false
% 
% 'PC_win'         nr of smaller time windows (time window of plot_win will
% be divided to PC_WIN nr of smaller windows to use for
% SPC calculation) 1 | 2
% 
% DS     spike or trial nr. downsampled
% 'no' | 'spike'
% 
% PLOT_WIN          1x2 vector, time window relative to event timestamp in
% sec, to use for SPC calculation
%
%
% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu
   


global cell_dir group_dir

% Chan label
if strcmp(rectype,'LFP')&& chanmean
    chtit = ['chanmean_' pr.subregion 'STN'];
elseif strcmp(rectype,'LFP')&& ~chanmean
    chtit = 'by-channel';
elseif strcmp(rectype,'EEG')
    chtit = 'F4';
end

% Freq label
if dominantfreq && ~contains(fr_name,'dom')
    frnm = ['dom_' fr_name];
else
    frnm = [fr_name];
end
frnm  = [frnm '_' num2str(plot_win)];
%% Load PC results (contains phase distribution parameters for each cell)

PC_results_dir = fullfile(cell_dir,'PC',[rectype '_phase_coupling'],chtit, event, event,frnm);
load(fullfile(PC_results_dir, ['PC_results_ds' ds '_' num2str(PC_win) 'win.mat']));

% Get PC  values
resvect = structfun(@(x)  x.ResVect, PC_results.Hilb_PC);
all_mrl = abs(resvect);
all_phas = angle(resvect);
all_cellids = fieldnames(PC_results.Hilb_PC );
for ic = 1:length(all_cellids)
    all_cellids{ic}(end-1) = '.';
end
    
% Only cells with significant PC

if contains(group,'resp')
    load(fullfile(group_dir,'RespCells.mat'));
    respevent =  group(1:strfind(group,'resp')-2);
    
  
    pdcellids{1} = RespCells.(respevent).none.Activ;
    pdcellids{2} = RespCells.(respevent).none.Inhib;
   
    
    resptyp = {'Activ','Inhib'};
    
elseif contains(group,'pred')
    load(fullfile(group_dir,'PredCells.mat'));
    respevent =   group(1:strfind(group,'pred')-2);
        
    pdcellids{1} = PredCells.([respevent '_StopPartition']).none.F_smaller_than_S;
    pdcellids{2} = PredCells.([respevent '_StopPartition' ]).none.F_bigger_than_S;
    
    resptyp = {'Fl_S','F_lS'};
elseif contains(group,'signPC');

    ray_sign = structfun(@(x) x.Ray_P<=0.05,PC_results.Hilb_PC);
 
    pdcellids{1} = all_cellids(ray_sign);
    respevent = 'sign PC'; resptyp = {};
    
elseif contains(upper(group),'UA')
    ray_sign = structfun(@(x) x.Ray_P<=0.05,PC_results.Hilb_PC);
    props = get_prop('SUA',all_cellids);
    
    cix1 = ray_sign+props==2;
    pdcellids{1} = all_cellids(cix1);
    cix2 = ray_sign+(~props)==2;
    pdcellids{2} = all_cellids(cix2);
      
    respevent = group; resptyp = {'SUA','MUA'};
    
else
    pdcellids{1}  = findcell;
    respevent = 'all'; resptyp = {};
end

[pdcellids{1},inx{1}] = intersect(all_cellids, pdcellids{1},'stable');
mrl{1} = all_mrl(inx{1});
phas{1} = all_phas(inx{1});

if length(pdcellids)>1
    [pdcellids{2},inx{2}] = intersect(all_cellids, pdcellids{2},'stable');
    mrl{2} = all_mrl(inx{2});
    phas{2} = all_phas(inx{2});
end
varargout{1} = respevent;
varargout{2} = resptyp;

end