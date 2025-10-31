function comp_vals = PC_RT_correlation(param,event, fr_name, chanmean, rectype, group, dominantfreq,PC_win, ds,plot_win,unit_group )
%PC_RT_CORRELATION      RT/SSDp0.5 correaltion with spike-phase coupling (SPC)
%   PC_RT_CORRELATION(param,event, fr_name, chanmean, rectype, group, dominantfreq,PC_win, ds,plot_win )
%       correalates PARAM values with spike-phase coupling measures:
%       preferred phase and MRL, calculated around EVENT in PC_WIN windows,
%       within PLOT_WIN time window. 
% 
% Input parameters:
%     PARAM             char. array, behav. measure to correlate with SPC measures
%       'RT' | 'SSD05'
% 
%     EVENT             char. array of event label
% 
%     FR_NAME           char. array, label of frequency band of interest
% 
%     CHANMEAN          if true, channel averaged LFP data is used, if false
%                       LFP data is plotted channel-by-channel
%                       (relevant only for intraop LFP data) 
% 
%     RECTYPE           recording type to use for analyses ('EEG' | 'LFP')
% 
%     GROUP             cell array of unit groups;
%        'all' | 'signPC' | [EVENT ' resp'] | [EVENT ' pred'], where EVENT
%        is a char. array of an event label
% 
%     DOMINANTFREQ      if true, dominant frequency within predefined frequency
%                       band limits (FREQS) are used for PC calculation, see FIND_DOMINANT_FREQ_BANDS
%                       true | false 
% 
%     PC_WIN            nr of smaller time windows (time window of plot_win will
%                       be divided to PC_WIN nr of smaller windows to use for
%                       SPC calculation) 
% 
%     DS                spike or trial nr. downsampled
%                    'no' | 'spike' 
% 
%     PLOT_WIN          1x2 vector, time window relative to event timestamp in
%                       sec, to use for SPC calculation
%
%     UNIT_GROUP        'SUA' | 'MUA'
% 
%  See also: PC_CELL_LEVEL, POLYPREDCICALL_MOD, LINCIRC_CORR2

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


global cell_dir figdir_pd


[PC_results_dir, mrl, phas, cellids,frnm] = get_phas(event, fr_name, chanmean, rectype,group, dominantfreq,PC_win, ds,plot_win );
if iscell(phas); phas = phas{1};end;
if iscell(cellids); cellids = cellids{1}; end;
if iscell(mrl); mrl = mrl{1}; end;

phas = mod(phas,2*pi);

if contains(upper(unit_group),'UA')
    
    props = get_prop('SUA',cellids);
    if strcmpi(unit_group,'SUA')
        cellids = cellids(props);
        mrl = mrl(props);
        phas = phas(props);
    elseif strcmpi(unit_group,'MUA')
        cellids = cellids(~props);
        mrl = mrl(~props);
        phas = phas(~props);
    end
end

if strcmp(param,'RT')
    % Get RT values
    rtval = nan(length(cellids),1);
    for ic = 1:length(cellids)
        act_cellid = cellids{ic}; act_cellid(end-1) = '.';
        TE = loadcb(act_cellid,'TrialEvents');
        RT = TE.KeyPress1 - TE.StimulusOn;
        [~,rmx] = rmoutliers(RT);
        RT(rmx) = NaN;
        rtval(ic) = nanmedian(RT);
    end
else
    rtval = nan(length(cellids),1);
    ssrt_dir = fullfile(figdir_pd,'Behav','SSRT');
    load(fullfile(ssrt_dir,'intraop_stimoff','SSRT_results.mat'));
    
    for ic = 1:length(cellids)
        act_cellid = cellids{ic}; act_cellid(end-1) = '.';
         [patnm, sessnm, chan, ~] = cellid2tags(act_cellid);
         rtval(ic) = SSRT_results.(patnm).(['S' sessnm]).(param)/1000;
    end
end

%% Linear regression
frnmtit = fr_name; frnmtit(strfind(fr_name,'_')) = '-';

% MRL
fig = figure;

[p, R] = polypredcicall_mod(mrl,rtval,0.95,'robust',0.01); hold on;
hold on;
scatter(mrl,rtval,[],'k','filled');
xlabel('MRL'); ylabel([param ' (s)']);
title({param ' vs MRL', ['Cells sign. coupled to ' frnmtit],[num2str(plot_win) ' sec around ' event]})

figdir =  fullfile(cell_dir,'PC',[param '_PC_correlation'],[rectype '_PC']); if ~isfolder(figdir); mkdir(figdir); end;
fnm = fullfile(figdir,['MRL_'  group '_' frnm '_' event '_ds' ds '_' num2str(PC_win) 'win_' num2str(plot_win)]);

if ~isempty(unit_group);
    fnm = [fnm '_' unit_group];
end
saveas(fig, [fnm '.fig'])
saveas(fig, [fnm '.jpg'])
saveas(fig, [fnm '.pdf'])
close(fig);

%% Fit mixture of gaussian distr.


 maxiter = 2000;
    Options.Display = 'off';
       Options.MaxIter = maxiter;
    Options.TolFun = 1e-6;
    compnr = 5;
%     colors = [0 0 0; 0 0 .5; 1 0.5 0; .5 .5 .5; 0 .5 0];
    
comp_vals = cell(compnr,compnr);


fig = figure;
rng(1);
for j = 1:compnr
    subplot(2,compnr,j)
GMModel = fitgmdist([mrl rtval],j,'start','plus','replicates',1,'options',Options);

idx = cluster(GMModel,[mrl rtval]);

Sc = nan(1,j);
for cln = 1:j
    Sc(cln) = scatter(mrl(idx==cln),rtval(idx==cln),'filled'); hold on;
    comp_vals{j,cln} = [mrl(idx==cln) rtval(idx==cln)];
end
xlabel('MRL'); ylabel([param ' (s)']);

legend(Sc, arrayfun(@(x) ['C.nr' num2str(x)],1:j,'UniformOutput',0),'Autoupdate','off');
setmyplot_balazs(gca)
xL = xlim; yL = ylim;
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
fc2 =fcontour(gmPDF,[xL yL]);
bic(j) = GMModel.BIC;
aic(j) = GMModel.AIC;
% text(xL(1),yL(2)*.9,['AIC: ' num2str(GMModel.AIC) ', BIC: ' num2str(GMModel.BIC)])    
end
subplot(2,compnr,compnr+[1:compnr])
plot(bic,'Color','g'); hold on; plot(aic,'Color','r'); legend({'BIC','AIC'})
xticks([1:3]); xlabel('Nr. of components'); ylabel('BIC/AIC');
suptitle({['Gaussian Mixture Model'], ...
    ['Cells sign. coupled to ' frnmtit ', ' num2str(plot_win) ' sec around ' event]});

setmyplot_balazs(gca)
fnm = [fnm '_GMM'];

saveas(fig, [fnm '.fig'])
saveas(fig, [fnm '.jpg'])
saveas(fig, [fnm '.pdf'])
close(fig);
%%
% Phase
fig = figure;


[Rval, pval] = lincirc_corr2(rtval,phas);
phas_phas = cat(1,phas,phas+2*pi);
medRT_RT = repmat(rtval,[2,1]);

scatter(phas_phas,medRT_RT,[],'k','filled');
tcks = [0 pi/2 pi 3*pi/2 2*pi 2*pi+pi/2 3*pi 2*pi+3*pi/2 4*pi];
xticks(tcks); xticklabels(arrayfun(@num2str,tcks,'UniformOutput',0));
xlabel('Phase angle'); ylabel([param ' (s)']);
xlim([0 4*pi])

if pval<0.05; col = 'r'; else; col = 'k'; end;
xL = xlim; yL = ylim;
text(xL(1),yL(2)*0.9,['p=' num2str(pval) ', R=' num2str(Rval)]);


xlabel('Resultant vector'); ylabel([param ' (s)']);
title({param ' vs RV', ['Cells sign. coupled to ' frnmtit],[num2str(plot_win) ' sec around ' event]})
setmyplot_balazs(gca);

figdir =  fullfile(cell_dir,'PC',[param '_PC_correlation'],[rectype '_PC']); if ~isfolder(figdir); mkdir(figdir); end;
fnm = fullfile(figdir,['RVphase_' group '_' frnm '_' event '_ds' ds '_' num2str(PC_win) 'win']);

if ~isempty(unit_group);
    fnm = [fnm '_' unit_group];
end
saveas(fig, [fnm '.fig'])
saveas(fig, [fnm '.jpg'])
saveas(fig, [fnm '.pdf'])
close(fig);
end



%--------------------------------------------------------------------------
function [PC_results_dir, mrl, phas, pdcellids,frnm,varargout] = get_phas(event, fr_name, chanmean, rectype,group, dominantfreq,PC_win, ds, plot_win)

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