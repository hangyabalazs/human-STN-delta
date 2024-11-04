function SSDp05_compare(groups_nm,whichSSRT,stat,conditions)
%SSDP05_COMPARE compares SDDp0.5 and SSRT values.
%   SSDP05_COMPARE(groups_nm,whichSSRT,stat,conditions) plots 
%   the estimate of stop signal delay at 50% probability of stopping
%   (SSDp0.5) or estimated  stop signal reaction time (SSRT)
%   for specified patient groups and compares them across conditions.
%   SSDp0.5/ SSRT values are saved in 'SSRT_results.mat' (see SSRTime_PD.m).
%   Box- or bar plots are drawn, grouped by task condition.
%   Plots are saved to  [ figdir '\Behav\'] directory under a separate subfolder
%   corresponding to WHICHSSRT.

%
%   Input parameters:
%       GROUPS_NM       cell array of patient group labels
%               {'all'}, includes all patients
%               {'tremor-dominant','akinetic-rigid','mixed'}, draws separate
%               plots for each clinical group
%               {'RTdecrease','RTincrease'}, - draws separate
%               plots for each patient group based on preop-postop RT change
% 
%       WHICHSSRT   which type of SSRT task parameter to compare
%               'ssd05' | 'ssrt'
%
%       STAT            statistical test to use; statistically significant differences
%                       are marked with lines above boxplots + asterisk (*<0.05; **< 0.01; ***<0.001)
%               'nonparam'     calculates median RT/ Performance, applies Mann-Whitney U test,  draws boxplots.
%               'param'     calculates mean RT/ Perf, applies unpaired t-test, draws bar plots.
%
%       CONDITIONS      Nx2 cell array of task conditions to compare
%               -first column corresponds to recording time, second column to
%               DBS stimulation condition
%               -ex.:{'preop','stimoff';'intraop','stimoff';'postop','stimoff';'postop','stimon'};

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu



set(0, 'DefaultFigureRenderer', 'painters');
grnr = length(groups_nm);
grvals = cell(grnr,1);
for pgi = 1:length(groups_nm)
    
    grnm = groups_nm{pgi};
    
    if ~contains(grnm,'all')
        pats = clinical_groups({ grnm });
        pats = pats{1};
    else
        pats = {'pd01','pd02','pd03','pd04','pd05','pd07','pd08','pd10','pd12','pd13','pd14','pd15','pd16'};
    end
    
    
    
        ssrtmat = cell(1,size(conditions,1));
    for sides = 1:2
        switch sides; case 1; side = 'left'; case 2; side = 'right'; end;
        %%
        ssrts = SSRT_plot2(pats,side,conditions,'whichSSRT',whichSSRT,...
            'stat',stat,'grnm',grnm);
        
        grvals{pgi,sides} = ssrts;
        
        for c = 1:size(conditions,1)
            ssrtmat{c}(:,sides) = ssrts(:,c);
        end
    end
    
    if contains(groups_nm{1},'all')
        RT_allgroupin(ssrtmat,conditions,pats,whichSSRT)
    end
    
end


%% Compare groups

if ~contains(grnm,'all')
    ssd05_groupcomp(grvals(1:2,:),groups_nm(1:2),whichSSRT, conditions,stat)
end
end



%--------------------------------------------------------------------------
function ssrts = SSRT_plot2(pats,side,conditions,varargin)
%SSRT_PLOT2 plots and compares ssdp0.5/ssrt values across task conditions
% SSRT_PLOT2(pats,side,conditions,...) compares
%      SSDp0.5 or SSRT values of patients specified in PATS and tested side 
%      specified in SIDE, across conditions specified in CONDITIONS.
%       Figures are saved to [ figdir_pd \ 'Behav' \  [whichSSRT '_compare']] directory.
%
%       Required inputs:
%           PATS     cell array of patients (some conditions for some
%                       patients might be missing)

%           SIDE     char. array, tested side to plot
%
%           CONDITIONS  cell array of conditions to compare; corresponds to
%                       SESSLIST and RT cell arrays
%
%       Optional input arguments (see RT_perf_compare)
%           'grnm'        character array, label of patient group included in PATS (default: 'all')               

%           'stat'        statistical test to use; statistically significant differences
%                         are marked with lines above boxplots + asterisk (*<0.05; **< 0.01; ***<0.001)
%               'nonparam'     applies Kruskal Wallis test with post-hoc Tukey-Kramer test,  draws boxplots.
%               'param'      applies One-way analysis of variance, draws bar plots.

%           'whichSSRT'   which type of SSRT task parameter to compare
%               'ssd05' | 'ssrt' (default: ssd05)
%
%       Output parameters:
%           SSRTS     m x n matrix where m is the number of patients
%                     included and n is the number of conditions to plot
% 
% See also: CALC_RT
% 
% Johanna Petra Szabó, 10.2024


prs = inputParser;
addRequired(prs,'pats',@iscell);
addRequired(prs,'side',@ischar);
addRequired(prs,'conditions',@iscell);
addParameter(prs,'grnm','all',@ischar);
addParameter(prs,'stat','nonparam',@(x) iscell(x)|ischar(x));
addParameter(prs,'whichSSRT','ssd05',@ischar);
parse(prs,pats,side,conditions,varargin{:});
pr = prs.Results;

%%
global figdir_pd

ssrt_dir = fullfile(figdir_pd,'Behav','SSRT');

condnr = length(pr.conditions);


ssrt_dir2 = fullfile(ssrt_dir,[pr.whichSSRT '_compare']); if ~isdir(ssrt_dir2); mkdir(ssrt_dir2); end;
lab = arrayfun(@(x) [conditions{x,1} '_' conditions{x,2}],1:condnr,'UniformOutput', 0);

patnr = length(pats);
ssrts = nan(length(pats),condnr);
for rt = 1:condnr
    rectime = conditions{rt,1};
    tag = conditions{rt,2};
    
    load(fullfile(ssrt_dir,[rectime '_' tag],'SSRT_results.mat'));
    
    for pp = 1:patnr
        patnm = pats{pp};
        if ~isfield(SSRT_results,patnm);
            fprintf('%s has no SSRT vals\n',patnm); continue;
        end;
        fnms = fieldnames(SSRT_results.(patnm));
        sesinx = cellfun(@(x) ismember(side(1),x),fnms);
        if ~any(sesinx)
            fprintf('%s %s has no SSRT vals\n',patnm,side); continue;
        end
        ssrts(pp,rt) = SSRT_results.(patnm).(fnms{sesinx}).(pr.whichSSRT)/1000;
        
        if strcmp(pr.whichSSRT,'ssd05') &&  ssrts(pp,rt)<0
            ssrts(pp,rt)  = NaN;
            fprintf('%s %s %s %s negative SSDp0.5 excluded\n',patnm,pr.side,rectime,tag);
        end
    end
    
    
end

if strcmp(pr.stat,'nonparam')
    [p, anovatab, stats] = kruskalwallis(ssrts,lab,'off');
    [c,m,h,gnames] = multcompare(stats,'display','off');
    
elseif strcmp(pr.stat,'param')
    [p, anovatag, stats] =  anova1(ssrts,lab,'off');
    [c,m,h,gnames] = multcompare(stats,'display','off');
end



fig = figure;

if strcmp(pr.stat,'nonparam')
    boxplot(ssrts,lab); hold on;
    set_my_boxplot(gca)
elseif strcmp(pr.stat,'param')
    bar(nanmean(ssrts,1)); hold on;
    er = errorbar(nanmean(ssrts,1),nanstd(ssrts,1)); er.Color = 'k'; er.LineStyle = 'none';
    legend(er,'STD');
    xticks(1:length(ssrts)); xticklabels(lab)
end

% Draw lines where data is not missing
% cmap = jet(patnr);
hold on;
% for k = 1:condnr-1
%     for pp = 1:patnr
%         Ln(pp) = line([k k+1],ssrts(pp,k:k+1),'Color',cmap(pp,:),'LineWidth',2);
%     end
% end
for kk = 1:condnr
    for pp = 1:patnr
        % scatter(kk+randi([-100 100],1)*0.0015,ssrts(pp,kk),[],cmap(pp,:),'filled'); hold on;
        scatter(kk+randi([-100 100],1)*0.0015,ssrts(pp,kk),[],'k','filled'); hold on;
    end
end
% legend(Ln,pats)


title([ pr.whichSSRT ' - ' side ' side'])

if strcmp(pr.whichSSRT,'SSD05p')
    ylabel( [pr.whichSSRT ' (%)'] );
else
    ylabel( [pr.whichSSRT ' (s)'] );
end

bpax = gca;
setmyplot_balazs(gca)

yL = ylim;
if p<.05; col = 'r'; else; col = 'k'; end;
text(0.5,yL(2)*.9,['p = ' num2str(p)],'Color',col);


            
for ic = 1:size(c,1)
    comps{ic} = [num2str(c(ic,1)) '-' num2str(c(ic,2))];
    
    if c(ic,end)<0.001
        sign(1,ic) = 3;
    elseif  c(ic,end)<0.01
        sign(1,ic) = 2;
    elseif  c(ic,end)<0.05
        sign(1,ic) = 1;
    else
        sign(1,ic) = 0;
    end
    if sign(1,ic)~=0
        boxplot_astx(bpax,c(ic,1:2),sign(1,ic))
        set_my_boxplot(gca)
    end
end


set(fig,'Position',get(0,'Screensize'));

fnm = fullfile(ssrt_dir2,[side '_' pr.grnm '_conds' num2str(condnr) '_' pr.stat]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])

close(fig)
end



%--------------------------------------------------------------------------
function ssd05_groupcomp(grvals,groups,whichSSRT, conditions,stat)
% SSD05_GROUPCOMP(grvals,groups,whichSSRT, conditions,stat) compares
% RT values across patient groups in every conditions specified.
%
% Input parameters:
%   GRVALS      m x 1 cell array, where m is the number of groups to
%                  compare; each cell contains n x 1 cell array, where n is the number
%                  of tested side to plot; each cell contains a p x c matrix,
%                  where p is the number of patients included, c is the
%                  number of conditions to plot
% 
%   GROUPS      cell array of group labels
% 
%   WHICHSSRT   which type of SSRT task parameter to compare
%               'ssd05' | 'ssrt'
% 
%   CONDITIONS  cell array of conditions to analyse;
% 
%   STAT        statistical test to use; statistically significant differences
%                  are marked with lines above boxplots + asterisk (*<0.05; **< 0.01; ***<0.001)
%               'nonparam'     applies Mann-Whitney U test,  draws boxplots.
%               'param'      applies Two-sample t-test, draws bar plots.

global figdir_pd

condnr = size(conditions,1);


sides = {'left','right'};
sdnr = length(sides)+1;

rt_figdir = fullfile(figdir_pd, 'Behav','SSRT',whichSSRT); if ~isdir(rt_figdir); mkdir(rt_figdir); end;

rt_figdir_grcomp = fullfile(rt_figdir,[groups{1} '_' groups{2}]);if ~isdir(rt_figdir_grcomp); mkdir(rt_figdir_grcomp); end;


fig = figure;
spnr = 1;
for condi = 1:condnr
    for sidi = 1:sdnr
        subplot(condnr,sdnr,spnr);
        
        switch sidi
            case {1,2}
                gr1 = grvals{1,sidi}(:,condi);
                gr2 = grvals{2,sidi}(:,condi);
                sidtit = sides{sidi};
            case 3
                gr11 = grvals{1,1}(:,condi);
                gr12 = grvals{1,2}(:,condi);
                gr21 = grvals{2,1}(:,condi);
                gr22 = grvals{2,2}(:,condi);
                gr1 = [gr11; gr12]; gr2 = [gr21; gr22];
                sidtit = 'bothsides';
        end
        maxpatnr = max([length(gr1) length(gr2)]);
        grs_oneside = cat(2,[gr1; nan(maxpatnr-length(gr1),1)], [gr2; nan(maxpatnr-length(gr2),1)]);
        
        
        if strcmp(stat,'nonparam')
            
            % Stat
            [p, h, stats] = ranksum(gr1,gr2);
            
            
            boxplot(grs_oneside/1000,groups); hold on;
            datpointsplot(grs_oneside/1000,1:2,{},[])
            
        elseif strcmp(stat,'param')
            
            
            xbar1 = nanmean(grs_oneside/1000,1);
            xstd1 = nanstd(grs_oneside/1000,1);
            
            
            % Stat
            [h,p,ci,stats] = ttest2(gr1,gr2);
            
            
            bar(xbar1); hold on;
            er = errorbar(xbar1,xstd1); er.Color = 'k';er.LineStyle = 'none';
            xticks(1:length(xbar1)); xticklabels(groups)
            
        end
        
        
        hold on;
        datpointsplot(grs_oneside,1:2,groups,[]);
        
        setmyplot_balazs(gca)
        bpax = gca;
        title([conditions{condi,1} ' ' conditions{condi,2} ', ' sidtit])
        ylabel(whichSSRT)
        
        
        if p<0.001
            sign = 3;
        elseif  p<0.01
            sign = 2;
        elseif  p<0.05
            sign = 1;
        else
            sign = 0;
        end
        if sign~=0
            boxplot_astx(bpax,[1 2],sign)
        end
        
        
        
        spnr = spnr+1;
    end
end

set(fig,'Position',get(0,'Screensize'));

fnm = fullfile(rt_figdir_grcomp,[groups{1} '_' groups{2} '_conds' num2str(condnr) '_' stat]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.eps'])

close(fig)
end

