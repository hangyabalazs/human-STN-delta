function PC_pie(event,plot_win,fr_name,downsamp,PC_win,alpha,rectype,chanmean,subregion,PCdir,groups)
%PC_PIE Pie chart of significantly phase-coupled units
%  PC_PIE(event,plot_win,fr_name,downsamp,PC_win,alpha,rectype,chanmean,subregion,PCdir)
%   draws a pie chart respresenting the percentage of significantly
%   phase-coupled units from all detected units, based on spike-phase
%   coupling measures calculated around EVENT, in PC_WIN windows within
%   PLOT_WIN sec time window.
%
% Input parameters:
%     EVENT             char. array of event label
%
%     PLOT_WIN          1x2 vector, time window relative to event timestamp in
%                       sec, to use for SPC calculation
%
%     FR_NAME           char. array, label of frequency band of interest
%
%     DOWNSAMP          spike or trial nr. downsampled
%                    'no' | 'spike'
%
%     PC_WIN            nr of smaller time windows (time window of plot_win will
%                       be divided to PC_WIN nr of smaller windows to use for
%                       SPC calculation)
%     ALPHA             numeric, alpha level, above this threshold units
%                       are considered sign. coupled)
%
%     RECTYPE           recording type to use for analyses ('EEG' | 'LFP')
%
%     CHANMEAN          if true, channel averaged LFP data is used, if false
%                       LFP data is plotted channel-by-channel
%                       (relevant only for intraop LFP data)
%     SUBREGION         character array or cell array, uses channel data
%                       derived from listed STN subregions
%                       (relevant only for intraop LFP data)
%          'all' | 'Motor' | 'Limbic' | 'Associative'
%     PCDIR             path to save results and figures
%
% See also: PC_CELL_LEVEL


% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


global group_dir cell_dir

if contains(rectype,'LFP')&& chanmean
    chtit = ['chanmean_' subregion 'STN'];
elseif contains(rectype,'LFP')&& ~chanmean
    chtit = 'by-channel';
elseif strcmp(rectype,'EEG')
    chtit = 'F4';
end



resdir = fullfile(PCdir,chtit,event,event,[fr_name '_' num2str(plot_win)]);
load(fullfile(resdir,['PC_results_ds' downsamp '_' num2str(PC_win) 'win.mat']));
figdir = fullfile(resdir, 'signPC_pie_charts');  if ~isdir(figdir); mkdir(figdir); end;

all_cellids0 = fieldnames(PC_results.Hilb_PC);
all_cellids = all_cellids0;
for j = 1:length(all_cellids0)
    all_cellids{j}(end-1) = '.';
end
if PC_win==1
    ray_sign = structfun(@(x) x.Ray_P<=alpha,PC_results.Hilb_PC);
elseif C_win==2
    ray_sign1 = structfun(@(x) x.Ray_P(1)<=alpha,PC_results.Hilb_PC);
    ray_sign2 = structfun(@(x) x.Ray_P(2)<=alpha,PC_results.Hilb_PC);
    ray_sign = ray_sign1|ray_sign2;
end
PC_results = [];
signcell = sum(ray_sign);
nosigncell = sum(~ray_sign);

sign_cellids = all_cellids(ray_sign);
nosign_cellids = all_cellids(~ray_sign);


allall = findcell;
nopc_calc_units = allall( ~ismember(allall, all_cellids) );
nopc_calc = length( nopc_calc_units); % no PC calculated

fig = figure;
if isempty(groups)
    
    pie([signcell nosigncell]);
    hold on;
    legend({['Sign PC, n = ' num2str(signcell)], ['Not sign. PC, n = ' num2str(nosigncell)]});
    
    fnm = fullfile(figdir,'signPC_pie');
elseif strcmpi(groups,'SUA')
    cellids0 = fieldnames(PC_results.Hilb_PC); cellids = cell(length(cellids0),1);
    for j = 1:length(cellids0)
        cellids{j} = cellids0{j}; cellids{j}(end-1) = '.';
    end
    props = get_prop('SUA',cellids);
    
    signcell_sua = sum(ray_sign+props==2);
    signcell_mua = sum(ray_sign+(~props)==2);
    nosigncell_sua = sum(~ray_sign+props==2);
    nosigncell_mua = sum(~ray_sign+(~props)==2);
    
    props2 = get_prop('SUA',nopc_calc_units);
    
    nopc_calc_sua = sum(props2);
    nopc_calc_mua = sum(~props2);
    
    pie([signcell_sua signcell_mua nosigncell_sua nosigncell_mua nopc_calc_sua nopc_calc_mua]);
    hold on;
    legend({['Sign PC SUA, n = ' num2str(signcell_sua)],['Sign PC MUA, n = ' num2str(signcell_mua)],...
        ['Not sign. PC SUA, n = ' num2str(nosigncell_sua)],  ['Not sign. PC MUA, n = ' num2str(nosigncell_mua)],...
        ['No PC calc SUA, n = ' num2str(nopc_calc_sua)], ['No PC calc MUA, n = ' num2str(nopc_calc_mua)]},...
        'Location','best');
    
    fnm = fullfile(figdir,'signPC_pie_suamua');
elseif contains(groups,'resp')
    
    
    load(fullfile(group_dir,'RespCells.mat'));
    r_event = groups(1:strfind(groups,'resp')-2);
    
    act_cellids = RespCells.(r_event).none.Activ;
    inh_cellids = RespCells.(r_event).none.Inhib;
    noresp_cellids = all_cellids( ~ismember(all_cellids, cat(2, act_cellids,inh_cellids )) );
    
    pc_act = length(intersect(act_cellids, sign_cellids));
    pc_inh = length(intersect(inh_cellids, sign_cellids));
    pc_noresp =  length(intersect(noresp_cellids, sign_cellids));
    nopc_act = length(intersect(act_cellids, nosign_cellids));
    nopc_inh = length(intersect(inh_cellids, nosign_cellids));
    nopc_noresp = length(intersect(noresp_cellids, nosign_cellids));
    
    
    pie([pc_act pc_inh pc_noresp nopc_act nopc_inh nopc_noresp nopc_calc]);
    hold on;
    legend({['Sign PC, ' r_event '-activ, n = ' num2str(pc_act)],...
        ['Sign PC, ' r_event '-inhib, n = ' num2str(pc_inh)],...
        ['Sign PC, not ' r_event '-resp, n = ' num2str(pc_noresp)],...
        ['Not sign PC, ' r_event '-activ, n = ' num2str(nopc_act)],...
        ['Not sign PC, ' r_event '-inhib, n = ' num2str(nopc_inh)],...
        ['Not sign PC, not ' r_event '-resp, n = ' num2str(nopc_noresp)],...
        ['No PC calc, n = ' num2str(nopc_calc)],...
        },'Location', 'bestoutside');
    
    fnm = fullfile(figdir,['signPC_pie_' r_event '_resp']);
    
    
elseif contains(groups,'signPC')&&contains(groups,'LFP')
    
    PCdir = fullfile(cell_dir,'PC',['LFP_phase_coupling']);
    
    resdir = fullfile(PCdir,'by-channel',event,event,[fr_name '_' num2str(plot_win)]);
    load(fullfile(resdir,['PC_results_ds' downsamp '_' num2str(PC_win) 'win.mat']));
    
    LFPall_cellids0 = fieldnames(PC_results.Hilb_PC);
    LFPall_cellids = LFPall_cellids0;
    for j = 1:length(LFPall_cellids0)
        LFPall_cellids{j}(end-1) = '.';
    end
    if PC_win==1
        LFPray_sign = structfun(@(x) x.Ray_P<=alpha,PC_results.Hilb_PC);
    elseif PC_win==2
        LFPray_sign1 = structfun(@(x) x.Ray_P(1)<=alpha,PC_results.Hilb_PC);
        LFPray_sign2 = structfun(@(x) x.Ray_P(2)<=alpha,PC_results.Hilb_PC);
        LFPray_sign = LFPray_sign1|LFPray_sign2;
    end
    LFPsign_cellids = LFPall_cellids(LFPray_sign);
    

    EEGpc_LFPpc_cellids = intersect(LFPsign_cellids, sign_cellids);
    EEGpc_LFPpc = length(EEGpc_LFPpc_cellids);
     only_LFPpc = length(LFPsign_cellids) - EEGpc_LFPpc;
     only_EEGpc = length(sign_cellids) - EEGpc_LFPpc;
    
    
    pie([EEGpc_LFPpc only_LFPpc only_EEGpc]);
    hold on;
    legend({['Sign PC both EEG & LFP, n = ' num2str(EEGpc_LFPpc)],...
        ['Sign. PC only to LFP, n = ' num2str(only_LFPpc)],...
        ['Sign. PC only to EEG, n = ' num2str(only_EEGpc)]},...
        'Location', 'bestoutside');
    
    fnm = fullfile(figdir,'signPC_EEG_LFP_pie');
    
    
    
elseif contains(groups,'signPC')&&contains(groups,'beta')
        
    Ln = strfind(groups,'_');
    fr_name2 = groups((Ln(1)+1):end);
    resdir = fullfile(PCdir,'by-channel',event,event,[fr_name2 '_' num2str(plot_win)]);
    load(fullfile(resdir,['PC_results_ds' downsamp '_' num2str(PC_win) 'win.mat']));
    
    BETAall_cellids0 = fieldnames(PC_results.Hilb_PC);
    BETAall_cellids = BETAall_cellids0;
    for j = 1:length(BETAall_cellids0)
        BETAall_cellids{j}(end-1) = '.';
    end
    if PC_win==1
        BETAray_sign = structfun(@(x) x.Ray_P<=alpha,PC_results.Hilb_PC);
    elseif PC_win==2
        BETAray_sign1 = structfun(@(x) x.Ray_P(1)<=alpha,PC_results.Hilb_PC);
        BETAray_sign2 = structfun(@(x) x.Ray_P(2)<=alpha,PC_results.Hilb_PC);
        BETAray_sign = BETAray_sign1|BETAray_sign2;
    end
    BETAsign_cellids = BETAall_cellids(BETAray_sign);
    

    DELTApc_BETApc_cellids = intersect(BETAsign_cellids, sign_cellids);
    DELTApc_BETApc = length(DELTApc_BETApc_cellids);
     only_BETApc = length(BETAsign_cellids) - DELTApc_BETApc;
     only_DELTApc = length(sign_cellids) - DELTApc_BETApc;
    
    
    pie([DELTApc_BETApc only_BETApc only_DELTApc]);
    hold on;
    legend({['Sign PC both delta & beta, n = ' num2str(DELTApc_BETApc)],...
        ['Sign. PC only to beta, n = ' num2str(only_BETApc)],...
        ['Sign. PC only to delta, n = ' num2str(only_DELTApc)]},...
        'Location', 'bestoutside');
    
    fnm = fullfile(figdir,'signPC_DELTA_BETA_pie');
    
    
    
end


saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);
