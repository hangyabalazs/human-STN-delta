function avg_psth_stat(event,partition,resdir,plot_win,parttags,patgr_nm,yL,varargin)
%AVG_PSTH_STAT Compares PSTHs of partitioned trials
% AVG_PSTH_STAT(event,partition,resdir,plot_win) selects units responsive
% to EVENT. Partitions epochs around EVENT based on PARTITION criteria.
% Calculates average PSTHs  of partitioned trials for each unit. Performs
% permutation statistics with cluster correction to compare partitions.
%
% Input parameters:
%     ALIGNEVENT	character array, label of behavioral event
%
%     PARTITION     character array, label of partition criteria
%
%     RESDIR        result directory to save plots
%
%     PLOT_WIN      time limits (sec) to perform statistics and generate plots
%
%     PATGR_NM          char. array, label of patient group, if '' all
%                       patients are used
%
% See also: PARTITION_TRIALS, RESPONSESORTER_PD, STD_STAT

% Johanna Petra Szabó, Balázs Hangya, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

if nargin<8
    excl_beforeEv = [];
else
    excl_beforeEv= varargin{1};
end

resptypes = {'Activ','Inhib'};
isfig = true;

alphas = 0.05;

set(0, 'DefaultFigureRenderer', 'painters');

[mylabels, mycolors] = makeColorsLabels(@defineLabelsColors_pd,...
    {[partition(2:end) '=' num2str(parttags(1)) ],[partition(2:end) '=' num2str(parttags(2))]});


switch event
    case 'StimulusOn'
        bwin = [-2.5 -1];
    case 'StopSignal'
        if isempty(excl_beforeEv)
            bwin = [-3 -1.5];   % baseline window for MW-test
        else
            bwin = [-1.1 -0.1];
        end
    case 'KeyPress1'
        bwin = [-3 -1.5];
    case 'Feedback'
        bwin = [1.5 3];
end

%%

if contains(upper(patgr_nm), 'RT') || contains(upper(patgr_nm), 'SLOWER')
    [patgroups, groups_nm] = clinical_groups({patgr_nm},'intraop','stimoff');
    
    % Get cellids - patient groups
    allcells = {};
    pats = patgroups{1};
    for j = 1:length(pats)
        patcells = findcell('rat',pats{j});
        allcells = cat(2,allcells,patcells);
    end
elseif contains(upper(patgr_nm), 'UA')
    cellids = findcell;
    props = get_prop('SUA',cellids);
    if strcmp(upper(patgr_nm),'SUA')
        allcells = cellids(props);
    elseif strcmp(upper(patgr_nm),'MUA')
        allcells = cellids(~props);
    end
    
else
    allcells = findcell;
end

time = -3:1/1000:3;
plottime = plot_win(1):1/1000:plot_win(2);
plotinx = dsearchn(time',plottime');


close all
cellids = cell(1,length(resptypes));
figure(1);
for rt = 1:length(resptypes)
    
    if contains(lower(patgr_nm), 'psth_stat')
        [isresp, cids] = get_prop('resp_psth_stat',allcells, 'Events2Align',{event});
    else
        [isresp, cids] = get_prop('resp',allcells, 'Events2Align',{event});
    end
    
    if strcmp(resptypes{rt},'Activ')
        cellids{rt} = cids(isresp==1);
    elseif  strcmp(resptypes{rt},'Inhib')
        cellids{rt} = cids(isresp==-1);
    end
    
    figure(2)
    psth_n = [];
    [psth_n] = norm_psth_map1(cellids{rt} ,resptypes{rt},event,...
        'baseline','indiv','basl_psth',{},'bwin',bwin, 'parts',partition,'cLim',[-6 6],...
        'bindex',[],'isfig',isfig,'sigma',0.08,'parttags',parttags,'excl_beforeEv',excl_beforeEv);
    
    
    % Downsample for stat
    
    xL = xlim;
    xplotL = interp1([-3 3],xL,plot_win);
    
    subplot(211)
    set(gca,'XLim',xplotL)
    subplot(212)
    set(gca,'XLim',xplotL)
    
    
    
    fnm =[partition(2:end) '_psth_MAP_' event '_' resptypes{rt} '_' patgr_nm];
    if ~isempty(excl_beforeEv)
        fnm = [fnm '_ExcBef_' excl_beforeEv];
    end
    saveas(figure(2),fullfile(resdir,[fnm '.jpg']))
    saveas(figure(2),fullfile(resdir,[fnm '.fig']))
    %     saveas(figure(2),fullfile(resdir,[fnm '.emf']))
    close(figure(2));
    
    
    
    
    figure(1);
    subplot(2,1,rt)
    [pL1,pL2,txt1,txt2] = avg_psth_plot(psth_n{1}(:,plotinx),psth_n{2}(:,plotinx), event, ...
        [mycolors{1}; mycolors{2}],plot_win);
    legend([pL1(end) pL2(end)],mylabels,'AutoUpdate','off','Location','best outside');
    if ~isempty(yL)
    ylim(yL);
    end
    
    statdat = {psth_n{1}';psth_n{2}'};
    
    [pcond, ~, ~,~,~,~] = std_stat(statdat,'condstats','on',...
        'mode','fieldtrip','fieldtripmethod','montecarlo','mcorrect','cluster',...
        'fieldtripnaccu',1000,...
        'fieldtripalpha',alphas,'fieldtripclusterparam',{'clusterstatistic','maxsize'});
    
    %
    %         [pcond, ~, ~,~,~,~] = std_stat(statdat,'condstats','on',...
    %             'mode','fieldtrip','fieldtripmethod','montecarlo','mcorrect','cluster',...
    %               'fieldtripnaccu',1000,'fieldtripalpha',alphas);
    setmyplot_balazs(gca);
    if isempty(pcond); continue; end;
    sign_p = pcond{1};
    
    hold on
    draw_signifpatch(time,sign_p,[0.5 0.5 0.5])
    
    
end

fnm =[partition(2:end) '_psth_AVG_' event '_resp_' patgr_nm];
if ~isempty(excl_beforeEv)
    fnm = [fnm '_ExcBef_' excl_beforeEv];
end
saveas(figure(1),fullfile(resdir,[fnm '.jpg']))
saveas(figure(1),fullfile(resdir,[fnm '.fig']))
saveas(figure(1),fullfile(resdir,[fnm '.pdf']))
close(figure(1));

end