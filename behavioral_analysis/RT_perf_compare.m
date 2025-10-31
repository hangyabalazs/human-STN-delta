function RT_perf_compare(groups_nm,conditions,stat,remoutliers, indivfig,avgfig,addlines,addpoints,rtplots,perfplots)
%RT_PERF_COMPARE    Comparison of reaction time (RT), performance and stop performance.
%   RT_perf_compare(groups_nm,conditions,...) calculates RT and/or performance
%   for specified patient groups and compares them across conditions.
%   Box- or bar plots are drawn, grouped by task condition.
%   Plots are saved to  [ figdir '\Behav\'] directory under a separate subfolder.
%   All_RTs.mat/All_Perfs.mat structs are saved to the same subdirectory.
%
%   Input parameters:
%       GROUPS_NM       cell array of patient group labels
%               {'all'}, includes all patients
%               {'tremor-dominant','akinetic-rigid','mixed'}, draws separate
%               plots for each clinical group
%               {'RTdecrease','RTincrease'}, - draws separate
%               plots for each patient group based on preop-postop RT change
%
%       CONDITIONS      1x2 cell array of task conditions to compare
%               -first column corresponds to recording time, second column to
%               DBS stimulation condition
%               -ex.:{'preop','stimoff';'intraop','stimoff';'postop','stimoff';'postop','stimon'};
%
%       STAT            statistical test to use; statistically significant differences
%                       are marked with lines above boxplots + asterisk (*<0.05; **< 0.01; ***<0.001)
%               'nonparam'     calculates median RT/ Performance, 
%                              applies Mann-Whitney U test/ Kruskal Wallis test with post-hoc Tukey Kramer test,
%                              draws boxplots.
%               'param'     calculates mean RT/ Perf, 
%                           applies Two sample t-test/ One way analysis of variance with post-hoc Tukey Kramer test,        
%                           draws bar plots.
%
%       REMOUTLIERS     specifies which outliers to remove
%                'indiv'    removes trials with outlier RT/Perf values for each patient
%                'avg'      removes patients with outlier RT/Perf values
%                'both'     both 'indiv' and 'avg'
%                'none'     does not remove any outliers

%       INDIVFIG         if true, performs comparison and draws plots for each individual patient
%               true | false
%
%
%       AVGFIG           if true, performs comparison and draws plots for
%                     pooled patient data, using patient averages
%               true | false

%       ADDLINES         if true, adds lines with distinct colour for each patient,
%                        connecting datapoing corresponding to the same patient
%                        across task conditions
%               true | false
%
%       ADDPOINTS        if true, adds scatter plots overlaid to box- or bar plots
%               true | false
%
%       RTPLOTS          calculates and compares RTs
%               true | false
%
%       PERFPLOTS        calculates and compares performance values
%               true | false
%

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu



global filesdir
 
%% Get RT values

condnr = size(conditions,1);
patgroups = {};
SessList = cell(condnr,1); RT = cell(condnr,1);

% Get patient groups if group label is specified
if ~contains(groups_nm{1},'all') 
    grnr = length(groups_nm);
    patgroups = clinical_groups(groups_nm);
else
    patgroups{1} = 'allpatients';
    grnr = 1;
end

grvals = cell(grnr,1);
for pgi = 1:grnr % loop over patient groups
    pats = patgroups{pgi};
    grnm = groups_nm{pgi};
    
    
    
    
    %% Load sess2analyse struct for each condition
    for coci = 1:condnr; % loop over conditions
        rectime = conditions{coci,1};
        condi = conditions{coci,2};
        
        SessList{coci} = getdata2analyse(filesdir, 'rectype','BEHAV',...
            'rectime',rectime,'patients', pats, 'side','bothside', 'condition',condi);
        
        TEnm = ['TrialEvents_nosync_' condi '.mat'];
        
        
        [RT{coci},Go_RT{coci},FAlarm_RT{coci},Hit{coci},NoStopTrials{coci},...
            StopTrials{coci},perf{coci},stopperf{coci}] = calc_RT(SessList{coci},TEnm,false);
        
    end
    
    allpats = unique({SessList{1}.patient});
    sides = unique({SessList{1}.side});
    
    
    
    %%% Compare RT across conditions
    if rtplots
        
        rtty = 'RT'; currRT = RT;
        % rtty = 'Go_RT'; currRT = Go_RT;
        % rtty = 'FAlarm_RT';  currRT = FAlarm_RT;
        
        if strcmp(stat,'param'); whichRT = [rtty 'mean']; else; whichRT = [rtty 'median']; end;
        
        
        % Stat. compare + plots of specified patient group
        allpat_cntr = RT_plots(allpats,grnm,sides,SessList,currRT,conditions,...
            'whichRT',whichRT,'stat',stat,'indivfig',indivfig,'avgfig',avgfig,...
            'remoutliers',remoutliers,'addlines',addlines,'addpoints',addpoints);
        
        grvals{pgi,1} = allpat_cntr;
        
        % Plot all patients marking which patient groups they are assign for
        if contains(groups_nm{1},'all')
            RT_allgroupin(allpat_cntr,conditions,allpats,whichRT)
        end
    end
    
    
    %%% Compare performance across conditions
    if perfplots
        
        % Stat. compare + plots of specified patient group
        % Performance 
        perf_sess = Perf_plots(allpats,grnm,sides,SessList,perf,conditions,'whichPerf','Perf','stat',stat,'indivfig',indivfig,'avgfig',avgfig,...
            'addlines',addlines,'addpoints',addpoints);
        
        % Stop performance
        stopperf_sess = Perf_plots(allpats,grnm,sides,SessList,stopperf,conditions,'whichPerf','StopPerf','stat',stat,'indivfig',indivfig,'avgfig',avgfig,...
            'addlines',addlines,'addpoints',addpoints);
        
        % Plot all patients marking which patient groups they are assign for
        if contains(groups_nm{1},'all')
            RT_allgroupin(perf_sess,conditions,allpats,'Perf')
            RT_allgroupin(stopperf_sess,conditions,allpats,'StopPerf')
        end
    end
    
    
end



if ~contains(groups_nm{1},'all')
    
%%% Compare RT across clinical groups        
    rtty = 'RT';
    % rtty = 'Go_RT'; currRT = Go_RT;
    % rtty = 'FAlarm_RT';  currRT = FAlarm_RT;
    
    
    if strcmp(stat,'param'); whichRT = [rtty 'mean']; else; whichRT = [rtty 'median']; end;
    
    if rtplots
        RT_compgroups(grvals,groups_nm(1:2),whichRT,conditions,stat)
    end
    
    
%%% Compare performance across groups
    if perfplots
        
        % Performance
        perf_comparegroups(groups_nm(1:2),conditions);
        
    end
end


end



%--------------------------------------------------------------------------
function allpat_cntr = RT_plots(allpats,grnm,sides,SessList,RT,conditions,varargin)
%RT_PLOTS compares and plots reaction time values across task conditions
% RT_plots(allpats,grnm,sides,SessList,RT,conditions,...) compares
%       reaction time (RT) values of patients specified in ALLPATS, across
%       conditions specified in CONDITIONS.
%       Figures are saved to [ figdir_pd '\Behav\RT'] directory.
%
%       Required inputs:
%           ALLPATS     cell array of patients (some conditions for some
%                       patients might be missing)
%
%           GRNM        character array, label of patient group included in ALLPATS
%
%           SIDES       cell array of task sides included in SESSLIST ('left' | 'right')
%
%           SESSLIST    cell array; each cell contains a struct with the list of data
%                       directories and other task information (patient code, side,
%                       stimulation condition, recording time and type for a given task condition
%                       (see RT_perf_compare, getdata2analyse)
% 
%           RT          cell array corresponding to task conditions, as in SESSLIST;
%                       each cell contains a cell array corresponding to
%                       the task settings listed in structs, stored in the
%                       cell arrays of SESSLIST; each cell contains the RT values for all
%                       trials corresponding to that task setting.
% 
%           CONDITIONS  cell array of conditions to compare; corresponds to
%                       SESSLIST and RT cell arrays
%
%       Optional input arguments (see RT_perf_compare)
%
%       Output parameters:
%           ALLPAT_CNTR    n x 1 cell array, where n is the number
%                          of task conditions to analyse; each cell contains a p x 2 matrix,
%                          where p is the number of patients included; the first column
%                          corresponds to left side task, second column to right side task
%                          results
% 
% See also: CALC_RT
% 
% Johanna Petra Szabó, 10.2024

prs = inputParser;
addRequired(prs,'allpats',@iscell);
addRequired(prs,'grnm',@ischar);
addRequired(prs,'sides',@iscell);
addRequired(prs,'SessList',@iscell);
addRequired(prs,'RT',@iscell);
addRequired(prs,'conditions',@iscell);
addParameter(prs,'whichRT','RT',@(x) iscell(x)|ischar(x));
addParameter(prs,'stat','nonparam',@(x) iscell(x)|ischar(x));
addParameter(prs,'addlines',false,@islogical);
addParameter(prs,'addpoints',true,@islogical);
addParameter(prs,'indivfig',true,@islogical);
addParameter(prs,'avgfig',true,@islogical);
addParameter(prs,'remoutliers','indiv',@ischar);

parse(prs,allpats,grnm,sides,SessList,RT,conditions,varargin{:});
pr = prs.Results;

global  figdir_pd rootdir
sides = unique({SessList{1}.side});
condnr = size(conditions,1);

lab = arrayfun(@(x) [conditions{x,1} ' ' conditions{x,2}],1:condnr,'UniformOutput', 0);

% load(fullfile(rootdir,'STN_loc','SubRegList.mat'));
% subreg_names = fieldnames(SubRegList);
%%

rt_figdir = fullfile(figdir_pd, 'Behav',pr.whichRT); if ~isdir(rt_figdir); mkdir(rt_figdir); end;



if pr.indivfig
    rt_figdir_gr = fullfile(rt_figdir,[grnm '_patfigs']);if ~isdir(rt_figdir_gr); mkdir(rt_figdir_gr); end;
end

rt_sess  = cell(condnr,1);
% sr_cods = {};
sign = [];
for j = 1:size(pr.allpats,2)
    patnm = pr.allpats{j};
    for jj = 1:2
        side = sides{jj};
        
        %% Align values for each patient {condition}{patient, side}
        %   (fill missing patients/side/condition data)
        for coci = 1:condnr
            
            rt_sess_pat{coci} = patside_pos(SessList{coci},patnm,side,RT{coci})';
            if strcmp(pr.remoutliers,'both') || strcmp(pr.remoutliers,'indiv');
                [~,rmx] = rmoutliers(rt_sess_pat{coci});
                rt_sess_pat{coci}(rmx) = NaN;
            end
            [rt_sess{coci}{j,jj}] =  rt_sess_pat{coci};
        end; clear coci;
        
        maxtrnr = max(cellfun(@length,rt_sess_pat));
        
        
        rt_sess_pat = cellfun(@(x) [x; nan(maxtrnr-length(x),1)],rt_sess_pat,'UniformOutput',0);
        
        %% Plot +  stat for individual patients
        if pr.indivfig
            kwx = cat(2,rt_sess_pat{:});
            
            kwx(find(kwx==0)) = NaN;
            gooddat = logical([condnr,1]);
            for jk = 1:condnr
                gooddat(jk) = any(~isnan(kwx(:,jk)));
            end
            
            
            
            % Stat
            try
                c = nan(condnr, 6); m = nan(condnr,2);
                if strcmp(pr.stat,'nonparam')
                    
                    if condnr>2
                        [pva,~,stats] = kruskalwallis(kwx(:,gooddat),lab(gooddat),'off');
                        [c,~,~,gnames] = multcompare(stats,'display','off');
                    else
                        [pva, ~, stats] = ranksum(kwx(:,1),kwx(:,2));
                        c(1,1) = 1; c(1,2) = 2; c(1,3) = pva;
                        
                    end
                    
                elseif strcmp(pr.stat,'param')
                    
                    [pva,~,stats] = anova1(kwx(:,gooddat),lab(gooddat),'off');
                    [c,~,~,~] = multcompare(stats,'display','off');
                    
                end
            catch
                fprintf('No stat could be computed %s %s\n',patnm,side)
                continue
            end
            
            
            % Boxplot
            fig = figure;
            if strcmp(pr.stat,'nonparam')
                boxplot(kwx(:,gooddat),'label',lab(gooddat)); hold on;
                delete(findobj(gcf,'tag','Outliers'))
                set_my_boxplot(gca)
                
                if pr.addpoints
                    datpointsplot(kwx,gooddat,lab,[0 0 0])
                end
                
            elseif strcmp(pr.stat,'param')
                std = cellfun(@nanstd,rt_sess_pat(gooddat));
                mean = cellfun(@nanmean,rt_sess_pat(gooddat));;
                
                
                bar(mean); hold on;
                er = errorbar(mean,std); er.Color = 'k';er.LineStyle = 'none';
                xticks(1:length(mean)); xticklabels(lab(gooddat))
                
            end
            %             ylim([0 3.5])
            bpax = gca;
            
            ylabel(pr.whichRT)
            title([patnm side])
            
            yL = ylim;
            if pva<.05; col = 'r'; else; col = 'k'; end;
            text(0.5,yL(2)*.9,['p=' num2str(pva)],'Color',col);
            
            for ic = 1:size(c,1)
                comps{ic} = [num2str(c(ic,1)) '-' num2str(c(ic,2))];
                comp_labs{ic} = {gnames{c(ic,1)}, gnames{c(ic,2)}};
                
                pvalc = c(ic,end);
                if pvalc<0.001
                    sign(j,jj,ic) = 3;
                elseif  pvalc<0.01
                    sign(j,jj,ic) = 2;
                elseif  pvalc<0.05
                    sign(j,jj,ic) = 1;
                else
                    sign(j,jj,ic) = 0;
                end
                
                
                if sign(j,jj,ic)~=0
                    boxplot_astx(bpax,c(ic,1:2),sign(j,jj,ic))
                end
                
                Pat_stat.(patnm).(side)(ic).Compared_groups =  comp_labs{ic};
                Pat_stat.(patnm).(side)(ic).Pval =  pvalc;
                Pat_stat.(patnm).(side)(ic).Diff_estimate =  c(ic,4);
            end
            Pat_stat.(patnm).(side)(ic).Stats =  stats;
            saveas(fig,fullfile(rt_figdir_gr,[patnm '_' side '.jpg']));
            saveas(fig,fullfile(rt_figdir_gr,[patnm '_' side '.pdf']));
            savefig(fig,fullfile(rt_figdir_gr,[patnm '_' side]));
            close(fig)
            
        end
    end
end
if pr.indivfig
    save(fullfile(rt_figdir_gr,'Pat_stat.mat'),'Pat_stat')
end

% %% Significance map
% ffig = figure;
% subplot(211)
% imagesc(squeeze(sign(:,1,:))); colorbar; colormap(winter)
% title([sides{1} ' side'])
% xticks(1:size(c,1))
% xticklabels(comps)
% ylabel('Patients')
%
% subplot(212)
% imagesc(squeeze(sign(:,2,:))); colorbar; colormap(winter)
% title([sides{2} ' side'])
% xticks(1:size(c,1))
% xticklabels(comps)
% ylabel('Patients')
% suptitle('Significance maps')
%
% saveas(ffig,fullfile(rt_figdir,'Sign_mapr.jpg'))
% saveas(ffig,fullfile(rt_figdir,'Sign_mapr.fig'))
% close(ffig)
%% Save aligned RT values

rti = 1;
for j = 1:size(pr.allpats,2)
    patnm = pr.allpats{j};
    for jj = 1:2
        side = sides{jj};
        
        All_RTs(rti).patient = patnm;
        All_RTs(rti).side = side;
        for coci = 1:condnr
            All_RTs(rti).([conditions{coci,1} '_' conditions{coci,2}]) = rt_sess{coci}{j,jj};
        end
        rti = rti+1;
        
    end
end

save(fullfile(rt_figdir,['All_RTs_' grnm '.mat']),'All_RTs');
%% Calculate patient mean/median


sdnr = length(sides);

% Calculate central values + remove outliers if necessary
for coci = 1:condnr
    if strcmp(pr.stat,'nonparam')
        allpat_cntr{coci} = cellfun(@nanmedian,rt_sess{coci});
        if strcmp(pr.remoutliers,'both') || strcmp(pr.remoutliers,'avg')
            %             [~,remix] = rmoutliers(allpat_cntr{coci});
            remix = allpat_cntr{coci} > 2.5;
            allpat_cntr{coci}(remix) = NaN;
        end
    elseif strcmp(pr.stat,'param')
        for coci = 1:condnr
            allpat_cntr{coci} = cellfun(@nanmean,rt_sess{coci});
            if strcmp(pr.remoutliers,'both') || strcmp(pr.remoutliers,'avg');
                % [~,remix] = rmoutliers(allpat_cntr{coci});
                remix = allpat_cntr{coci} > 2.5;
                allpat_cntr{coci}(remix) = NaN;
            end
        end
    end
    
    allpat_cntr{coci}(allpat_cntr{coci}==0) = NaN;
end


%% Stat + plots for each side (using patient-averages)
if pr.avgfig
    fig = figure;
    for sidi = 1:sdnr
        
        x_oneside = [];
        for coci = 1:condnr
            x_oneside = [x_oneside allpat_cntr{coci}(:,sidi)];
        end
        
        
        subplot(sdnr,1,sidi);
        
        
        if strcmp(pr.stat,'nonparam')
            
            % Stat
            if condnr>2
                [p1,anovatab1,stats1] = kruskalwallis(x_oneside,lab,'off');
                [c1,m1,h1,gnames1] = multcompare(stats1,'display','off');
            else
                [p1, ~, stats1] = ranksum(x_oneside(:,1),x_oneside(:,2));
                c1(1,1) = 1; c1(1,2) = 2; c1(1,3) = p1;
            end
            
            
            
            boxplot(x_oneside,lab); hold on;
            set_my_boxplot(gca)
            
            if pr.addpoints
                datpointsplot(x_oneside,1:size(x_oneside,2),lab,[0 0 0])
            end
            
        elseif strcmp(pr.stat,'param')
            
            
            xbar1 = nanmean(x_oneside,1);
            xstd1 = nanstd(x_oneside,1);
            
            
            % Stat
            [p1,anovatab1,stats1] = anova1(x_oneside,lab,'off');
            [c1,m,h1,gnames1] = multcompare(stats1,'display','off');
            
            
            bar(xbar1); hold on;
            er = errorbar(xbar1,xstd1); er.Color = 'k';er.LineStyle = 'none';
            xticks(1:length(xbar1)); xticklabels(lab)
            
        elseif strcmp(pr.stat,'ranova')
            
            matx = mat2cell(x_oneside,13, [1 1 1 1]);
            anovlab = arrayfun(@(x)  ['cond' num2str(x)] ,[1:condnr],'UniformOutput',0);
            tab = table(allpats',matx{:},'VariableNames',{'Patients',anovlab{:}});
            WDis = table([1 2 3 4]','VariableNames',{'Cond'});
            
            rm = fitrm(tab,  [anovlab{1} '-' anovlab{end} '~Patients'], 'WithinDesign', WDis);
            ranovatbl = ranova(rm)
        end
        %         ylim([0 3.5])
        setmyplot_balazs(gca)
        bpax = gca;
        title(['RT across patients - ' sides{sidi} ' side'])
        ylabel(pr.whichRT)
        
        yL = ylim;
        if p1<.05; col = 'r'; else; col = 'k'; end;
        text(0.5,yL(2)*.9,['p=' num2str(p1)],'Color',col);
        
        for ic = 1:size(c1,1)
            comps{ic} = [num2str(c1(ic,1)) '-' num2str(c1(ic,2))];
            
            if c1(ic,end)<0.001
                sign(1,ic) = 3;
            elseif  c1(ic,end)<0.01
                sign(1,ic) = 2;
            elseif  c1(ic,end)<0.05
                sign(1,ic) = 1;
            else
                sign(1,ic) = 0;
            end
            if sign(1,ic)~=0
                boxplot_astx(bpax,c1(ic,1:2),sign(1,ic))
            end
        end
        
        
        
        if pr.addlines
            cmap = colormap(jet);
            cmap = cmap(1:floor(size(cmap,1)/length(pr.allpats)):end,:);
            for pati = 1:size(x_oneside,1);
                for compi = 1:size(x_oneside,2)-1
                    a = compi; b = compi+1;
                    Ln(pati) = line([a b],[x_oneside(pati,a) x_oneside(pati,b)],'Color',cmap(pati,:),'LineWidth',2);
                end
            end
            legend(Ln,pr.allpats)
        end
        
        
        
    end
    %     set(fig,'Position',get(0,'Screensize'));
    
    fnm = fullfile(rt_figdir,[grnm '_conds' num2str(condnr) '_lines_' num2str(pr.addlines) '_' pr.remoutliers 'rmout']);
    saveas(fig,[fnm '.jpg'])
    saveas(fig,[fnm '.fig'])
    saveas(fig,[fnm '.pdf'])
    
    close(fig)
    
end

end


%-------------------------------------------------------------------------------------
function perf_sess = Perf_plots(allpats,grnm,sides,SessList,perf,conditions,varargin)

%PERF_PLOTS compares and plots performance values across task conditions
% PERF_PLOTS(allpats,grnm,sides,SessList,perf,conditions,...) compares
%       reaction time (RT) values of patients specified in ALLPATS, across
%       conditions specified in CONDITIONS.
%       Figures are saved to [ figdir_pd '\Behav\Perf'] directory.
%
%       Required inputs:
%           ALLPATS     cell array of patients (some conditions for some
%                       patients might be missing)
%
%           GRNM        character array, label of patient group included in ALLPATS
%
%           SESSLIST    cell array; each cell contains a struct with the list of data
%                       directories and other task information (patient code, side,
%                       stimulation condition, recording time and type for a given task condition
%                       (see RT_perf_compare, getdata2analyse)
%           RT          cell array corresponding to task conditions, as in SESSLIST;
%                       each cell contains a cell array corresponding to
%                       the task settings listed in structs, stored in the
%                       cell arrays of SESSLIST; each cell contains the RT values for all
%                       trials corresponding to that task setting.
%           CONDITIONS  cell array of conditions to compare; corresponds to
%                       SESSLIST and RT cell arrays
% 
% See also: CALC_RT
%
% Johanna Petra Szabó, 10.2024


prs = inputParser;
addRequired(prs,'allpats',@iscell);
addRequired(prs,'grnm',@ischar);
addRequired(prs,'sides',@iscell);
addRequired(prs,'SessList',@iscell);
addRequired(prs,'perf',@iscell);
addRequired(prs,'conditions',@iscell);
addParameter(prs,'stat','nonparam',@(x) iscell(x)|ischar(x));
addParameter(prs,'whichPerf','Perf',@(x) iscell(x)|ischar(x));
addParameter(prs,'addlines',false,@islogical);
addParameter(prs,'addpoints',true,@islogical);
addParameter(prs,'indivfig',true,@islogical);
addParameter(prs,'avgfig',true,@islogical);

parse(prs,allpats,grnm,sides,SessList,perf,conditions,varargin{:});
pr = prs.Results;

global  figdir_pd
sides = unique({SessList{1}.side});
sdnr = length(sides);
condnr = size(conditions,1);

lab = arrayfun(@(x) [conditions{x,1} ' ' conditions{x,2}],1:condnr,'UniformOutput', 0);

% load(fullfile(rootdir,'STN_loc','SubRegList.mat'));
% subreg_names = fieldnames(SubRegList);
%%

perf_figdir = fullfile(figdir_pd, 'Behav',pr.whichPerf); if ~isdir(perf_figdir); mkdir(perf_figdir); end;


perf_sess  = cell(condnr,1);
sign = [];
for j = 1:size(pr.allpats,2)
    patnm = pr.allpats{j};
    for jj = 1:2
        side = sides{jj};
        
        %% Align values for each patient {condition}{patient, side}
        %   (fill missing patients/side/condition data)
        for coci = 1:condnr
            
            perf_sess{coci}(j,jj) = patside_pos(SessList{coci},patnm,side,perf{coci})';
            if perf_sess{coci}(j,jj)== 0;
                perf_sess{coci}(j,jj) = NaN;
            end
        end;
        
        
    end
end


%% Stat + plots for each side (using patient-averages)

fig = figure;
for sidi = 1:sdnr
    
    x_oneside = [];
    for coci = 1:condnr
        x_oneside = [x_oneside perf_sess{coci}(:,sidi)];
    end
    
    
    subplot(sdnr,1,sidi);
    
    
    if strcmp(pr.stat,'nonparam')
        
        % Stat
        if condnr>2
            [p1,anovatab1,stats1] = kruskalwallis(x_oneside,lab,'off');
            [c1,m1,h1,gnames1] = multcompare(stats1,'display','off');
        else
            [p1, h1, stats1] = ranksum(x_oneside(:,1),x_oneside(:,2));
            c1(1,1) = 1; c1(1,2) = 2; c1(1,3) = p1;
            
        end
        
        
        boxplot(x_oneside,lab); hold on;
        set_my_boxplot(gca)
        
        if pr.addpoints
            datpointsplot(x_oneside,1:size(x_oneside,2),lab,[0 0 0])
        end
        
    elseif strcmp(pr.stat,'param')
        
        
        xbar1 = nanmean(x_oneside,1);
        xstd1 = nanstd(x_oneside,1);
        
        
        % Stat
        [p1,anovatab1,stats1] = anova1(x_oneside,lab,'off');
        [c1,m,h1,gnames1] = multcompare(stats1,'display','off');
        
        
        bar(xbar1); hold on;
        er = errorbar(xbar1,xstd1); er.Color = 'k';er.LineStyle = 'none';
        xticks(1:length(xbar1)); xticklabels(lab)
        
    end
    setmyplot_balazs(gca)
    bpax = gca;
    title([sides{sidi} ' side'])
    ylabel([pr.whichPerf ' (%)'])
    
    
    yL = ylim;
    if p1<.05; col = 'r'; else; col = 'k'; end;
    text(0.5,yL(2)*.9,['p=' num2str(p1)],'Color',col);
    for ic = 1:size(c1,1)
        comps{ic} = [num2str(c1(ic,1)) '-' num2str(c1(ic,2))];
        
        if c1(ic,end)<0.001
            sign(1,ic) = 3;
        elseif  c1(ic,end)<0.01
            sign(1,ic) = 2;
        elseif  c1(ic,end)<0.05
            sign(1,ic) = 1;
        else
            sign(1,ic) = 0;
        end
        if sign(1,ic)~=0
            boxplot_astx(bpax,c1(ic,1:2),sign(1,ic))
        end
    end
    
    %     ylim([30 100]);
    
    
    if pr.addlines
        cmap = colormap(jet);
        cmap = cmap(1:floor(size(cmap,1)/length(pr.allpats)):end,:);
        for pati = 1:size(x_oneside,1);
            for compi = 1:size(x_oneside,2)-1
                a = compi; b = compi+1;
                Ln(pati) = line([a b],[x_oneside(pati,a) x_oneside(pati,b)],'Color',cmap(pati,:),'LineWidth',2);
            end
        end
        legend(Ln,pr.allpats)
    end
    
    
    
end
% set(fig,'Position',get(0,'Screensize'));

fnm = fullfile(perf_figdir,[grnm '_conds' num2str(condnr) '_lines_' num2str(pr.addlines)]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])

close(fig)

% Cal percentage of perf<60;
allvals0 = cat(1,perf_sess{:});
allvals = allvals0(:);
lowperf_perf = ( sum(allvals>=60)*100 )/ length(allvals);
fprintf('%d perc of tests >=60 perc performance\n',lowperf_perf)


%%

rti = 1;
for j = 1:size(pr.allpats,2)
    patnm = pr.allpats{j};
    for jj = 1:2
        side = sides{jj};
        
        All_Perfs(rti).patient = patnm;
        All_Perfs(rti).side = side;
        for coci = 1:condnr
            All_Perfs(rti).([conditions{coci,1} '_' conditions{coci,2}]) = perf_sess{coci}(j,jj);
        end
        rti = rti+1;
        
    end
end

save(fullfile(perf_figdir,['All_Perfs_' grnm '.mat']),'All_Perfs');

end





% ----------------------------------------------------------------------------
function RT_compgroups(grvals,groups,whichRT, conditions,stat)
%RT_COMPGROUPS compares RT values across patient groups
% RT_compgroups(grvals,groups,whichRT, conditions,stat) compares
% RT values across patient groups in every conditions specified.
%
% Input parameters:
%       GRVALS     m x 1 cell array, where m is the number of groups to
%                  compare; each cell contains n x 1 cell array, where n is the number
%                  of task conditions to analyse; each cell contains a p x 2 matrix,
%                  where p is the number of patients included; the first column
%                  corresponds to left side task, second column to right side task
%                  results
%
%       GROUPS     cell array of group labels
%
%       WHICHRT    which type of RT to compare
%           'RT' | 'Go_RT' | 'FAlarm_RT' (see calc_RT.m)
%
%       CONDITIONS  cell array of conditions to analyse;
%
%       STAT       statistical test to use; statistically significant differences
%                  are marked with lines above boxplots + asterisk (*<0.05; **< 0.01; ***<0.001)
%               'nonparam'     calculates median RT/ Performance, applies Mann-Whitney U test,  draws boxplots.
%               'param'     calculates mean RT/ Perf, applies unpaired t-test, draws bar plots.
%
% Johanna Petra Szabó, 10.2024

global figdir_pd

condnr = size(conditions,1);


sides = {'left','right'};
% sdnr = length(sides)+1;
sdnr = length(sides);
rt_figdir = fullfile(figdir_pd, 'Behav',whichRT); if ~isdir(rt_figdir); mkdir(rt_figdir); end;

rt_figdir_grcomp = fullfile(rt_figdir,[groups{1} '_' groups{2}]);if ~isdir(rt_figdir_grcomp); mkdir(rt_figdir_grcomp); end;


fig = figure;
spnr = 1;
for condi = 1:condnr
    for sidi = 1:sdnr
        subplot(condnr,sdnr,spnr);
        %         subplot(sdnr,condnr,spnr);
        
        switch sidi
            case {1,2}
                gr1 = grvals{1}{condi}(:,sidi);
                gr2 = grvals{2}{condi}(:,sidi);
                sidtit = sides{sidi};
            case 3
                
                gr1 = [grvals{1}{condi}(:,1); grvals{1}{condi}(:,2)];
                gr2 = [grvals{2}{condi}(:,1); grvals{2}{condi}(:,2)];
                sidtit = 'bothsides';
        end
        maxpatnr = max([length(gr1) length(gr2)]);
        grs_oneside = cat(2,[gr1; nan(maxpatnr-length(gr1),1)], [gr2; nan(maxpatnr-length(gr2),1)]);
        
        
        if strcmp(stat,'nonparam')
            
            % Stat
            [p, h, stats] = ranksum(gr1,gr2);
            
            
            boxplot(grs_oneside,groups);
            set_my_boxplot(gca)
            
        elseif strcmp(stat,'param')
            
            
            xbar1 = nanmean(grs_oneside,1);
            xstd1 = nanstd(grs_oneside,1);
            
            
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
        ylabel([whichRT '(s)'])
        
        
        
        yL = ylim;
        if p<.05; col = 'r'; else; col = 'k'; end;
        text(0.5,yL(2)*.9,['p=' num2str(p)],'Color',col);
        
        
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

close(fig)

end



%--------------------------------------------------------------------------
function perf_comparegroups(groups_nm,conditions)
%PERF_COMPAREGROUPS compares performance values across patient groups
% PERF_COMPAREGROUPS(groups_nm,conditions) compares
% performance values across patient groups in every conditions specified.
%
% Input parameters:
%       GROUPS_NM    cell array of group labels
%
%       CONDITIONS  cell array of conditions to analyse;
%
% Johanna Petra Szabó, 10.2024


global filesdir figdir_pd

sides = {'left','right'};


for condi = 1:size(conditions,1)
    
    %%
    rectime = conditions{condi,1};
    tag = conditions{condi,2};
    
    patgroups = clinical_groups(groups_nm);
    grnr = length(groups_nm);
    
    perf = cell(1,grnr);
    stopperf = cell(1,grnr);
    
    for si = 1:length(sides)
        for pgi = 1:grnr
            pats = patgroups{pgi};
            grnm = groups_nm{pgi};
            
            
            
            
            %% Load sess2analyse struct for each condition
            
            
            SessList= getdata2analyse(filesdir, 'rectype','BEHAV',...
                'rectime',rectime,'patients', pats, 'side',sides{si}, 'condition',tag);
            
            TEnm = ['TrialEvents_nosync_' tag '.mat'];
            
            
            [~,~,~,~,~,~,perf{pgi,si},stopperf{pgi,si}] = calc_RT(SessList,TEnm,false);
            
        end
    end
    
    for k = 1:2 % 1 - perf, 2 - stop perf
        switch k;
            case 1;
                perfx = perf;
                whichperf = 'Perf';
                yL = [30 100];
            case 2;
                perfx = stopperf;
                whichperf = 'StopPerf';
                yL = [0 100];
        end;
        resdir = fullfile(figdir_pd,'Behav',whichperf,[groups_nm{1} '_' groups_nm{2}]);
        if ~isfolder(resdir); mkdir(resdir); end;
        
        perfx = cellfun(@cell2mat,perfx,'UniformOutput',0);
        pL = cellfun(@length,perfx);
        mL = max(pL(:));
        
        perfnns = cellfun(@(x) [x nan(1,mL-length(x))],perfx,'UniformOutput',0);
        
        fig = figure;
        %         for si = 1:3
        for si = 1:2
            %             subplot(3,1,si)
            subplot(2,1,si)
            switch si
                case 3
                    x1 = cat(1,perfnns{:,1})';
                    x2 = cat(1,perfnns{:,2})';
                    X = [x1;x2];
                    tit = 'bothsides'
                otherwise
                    X = cat(1,perfnns{:,si})';
                    tit = sides{si};
            end
            boxplot(X,groups_nm); hold on;
            set_my_boxplot(gca)
            
            ylim(yL)
            datpointsplot(X,1:grnr,{},[]);
            title(tit);
            ylabel([whichperf ' (%)'])
            %         [h,pval] = ttest2(X(:,1),X(:,2));
            
            [pval,h] = ranksum(X(:,1),X(:,2));
            
            if pval<0.05;col = 'r'; else; col = 'k'; end;
            yL = ylim;
            text(0.5,yL(2)*0.9,['p = ' num2str(pval)],'Color',col);
            
            setmyplot_balazs(gca)
        end
        %         set(gcf,'Position',get(0,'ScreenSize'));
        saveas(fig,fullfile(resdir,[rectime '_' tag '.jpg']))
        saveas(fig,fullfile(resdir,[rectime '_' tag '.fig']))
        saveas(fig,fullfile(resdir,[rectime '_' tag '.pdf']))
        close(fig);
    end
end
end

%--------------------------------------------------------------------------
function [rt_pat] = patside_pos(allsessions,patnm,side,RT)
%PATSIDE_POS assigns RT values to patient sessions
% [rt_pat] = patside_pos(allsessions,patnm,side,RT) finds RT value matching
% inputted patient (PATNM) & side (SIDE). If there is no such value, result = 0.

% Input parameters:
%   ALLSESSIONS   session information (see getdata2analyse)
%   PATNM         patient code
%   SIDE          side of experiment
%   RT            RT values correponding to one experimental condition
%
% Output parameter:
%   RT_PAT        RT values corresponding to PATNM, SIDE





patpos = find(ismember({allsessions.patient},patnm));
sesspos = find(ismember({allsessions(patpos).side},side));

if ~isempty(sesspos)&& patpos(sesspos)<=length(RT)
    rt_pat = RT{patpos(sesspos)};
else
    rt_pat = 0;
end
if isempty(rt_pat)
    rt_pat = 0;
end
end
