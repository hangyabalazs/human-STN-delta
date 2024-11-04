function unit_bursting(EventTypes)
%UNIT_BURSTING      Analysis of bursting properties of units
%   UNIT_BURSTING(EventTypes) 
%       -generates an auto-correlogram map (with units
%       sorted based on bursting indeces (BI) ), averaged auto-correlograms and
%       a pie chart of bursting and non-bursting units. An empiric threshold of
%       0.35 BI is used. 
%       -correlates BI with UPDRS scores and behavioral parameters (RT,
%       SSDp0.5)
%       -separates units to low- and high bursting groups based on median BI
%       -generates average PSTH plots (with standard error) for low- and
%       high bursting groups including all epochs around a selected event/
%       including epochs partitioned based on outcome of stopping.
%       The difference between Successful and Failed stop trials is tested
%       with permutation test with cluster based correction in both groups.
%       
%
% See also: AUTOCORR_PD, POLYPREDCICALL_MOD, STD_STAT, NORM_PSTH_MAP1, AVG_PSTH_PLOT

% Johanna Petra Szabó, Balázs Hangya, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

global figdir_pd
% Parameters
burst_color = rgb('maroon');
noburst_color = rgb('gray');
colors = [burst_color; noburst_color];

BIthr = 0.35; % threshold for bursting index
plot_win = [-0.05 0.05];

resdir = fullfile(figdir_pd,'Intraop_SP','Bursting'); if ~isdir(resdir); mkdir(resdir); end;
respcells = cellfun(@(x) cat(2,x,' resp'),EventTypes,'UniformOutput',0);




% Auto-correlogram map, sorted based on bursting index of units
ACGmap_sorted(BIthr,plot_win)


% Average auto-correlograms of bursting/ non-bursting units
avgACG_BI_nonBI(BIthr,plot_win,resdir,burst_color,noburst_color)


% Percentage of bursting units (pie chart)
bursting_pies(BIthr,colors,resdir)


% Correlation of bursting index with UPDRS scores
figdir = fullfile(figdir_pd,'Intraop_SP','Bursting','UPDRS_corrs'); if ~isfolder(figdir); mkdir(figdir); end;
updrs_bursting_corr(false, figdir)


% Correlation of bursting index with behavioral parameterss
for k = 1:2
    switch k; case 1; param = 'RT'; case 2; param = 'ssd05'; end;
    
    resdir_behav = fullfile(figdir_pd,'Intraop_SP','Bursting','behav_corrs'); if ~isdir(resdir_behav); mkdir(resdir_behav); end;
    
    bursting_RT(resdir_behav,param,false,'bothside')
end

% Peri-event time histogram (PSTH) of units separated by median bursting index

for ei = 1:length(EventTypes)
    
    respev = EventTypes(ei);
    
    burst_mediansplit('median',respev,resdir,'all',[.05, .01]) % PSTHs for all epochs around RESPEV
    burst_mediansplit('median',respev,resdir,'#StopPartition',[.05, .01]); % PSTHs for stop trials around RESPEV, partitioned based on outcome of stopping + statistics
end



end



%--------------------------------------------------------------------------
function bursting_pies(thr,colors,resdir)

pdcells = findcell;
brst = get_prop('bursting',pdcells);

fig = figure;

bc = brst>=thr; nbc = brst<thr;

pie([sum(bc) sum(nbc) sum(isnan(brst))]);
ax = gca;
ax.Colormap = [colors; 0 0 0];
legend({['BI>=' num2str(thr)],['BI<' num2str(thr)],'Not enough spikes'})

saveas(fig,fullfile(resdir,['pie_' num2str(thr) '.jpg' ]))
saveas(fig,fullfile(resdir,['pie_' num2str(thr) '.fig' ]))
saveas(fig,fullfile(resdir,['pie_' num2str(thr) '.pdf' ]))
close(fig);


end


%--------------------------------------------------------------------------
function ACGmap_sorted(BIthr,plot_win)

global cell_dir


% Time
timvec = -3:1/2000:3; timvec = timvec(1:end-1);
pltim = plot_win(1):1/2000:plot_win(2);
plotlims = dsearchn(timvec',plot_win');


% sort ACGs
load(fullfile(cell_dir,'ACG','high_delta','ACG_matrices_')); % burst index is the same for each freq band
[sortBI,sortix] = sort(BurstIndex,'descend');
sortix(isnan(sortBI)) = [];
sortBI(isnan(sortBI)) = [];
BIlim = find(sortBI<BIthr,1);

CCR_z = zscore(CCR')';
sortCCR = CCR_z(sortix,:);

% Plot
fig = figure;
pcolor(pltim,1:length(sortix),sortCCR(:,plotlims(1): plotlims(end))); colorbar;
shading interp; set(gca,'YDir','reverse');
% caxis(prctile( CCR_z(:), [1 99] ));
caxis([-2 8]);

colormap(pink)
hold on;
xL = xlim; yL = ylim;
line(xL,[BIlim BIlim],'Color','r','LineWidth',2);


cut_win = [-0.01 0.01];
line([cut_win; cut_win],repmat(yL',[1 2]),'Color','k','LineWidth',2)


xlabel('Time lag (s)');
ylabel('Units')

fnm = fullfile(cell_dir,'Bursting','ACG_map',['ACG_sorted_map' num2str(plot_win(1)) '_' num2str(plot_win(2))] );
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);
end



%--------------------------------------------------------------------------
function avgACG_BI_nonBI(BIthr,plot_win,resdir,burst_color,noburst_color)

global cell_dir



timvec = -3:1/2000:3; timvec = timvec(1:end-1);
pltim = plot_win(1):1/2000:plot_win(2);
plotlims = dsearchn(timvec',plot_win');

load(fullfile(cell_dir,'ACG','high_delta','ACG_matrices_')); % burst index is the same for each freq band
BI = BurstIndex>=BIthr;
nonBI = BurstIndex<BIthr;

BIdat = CCR(BI,plotlims(1):plotlims(end));
nonBIdat = CCR(nonBI,plotlims(1):plotlims(end));

fig = figure;
p1 = errorshade(pltim,nanmean(BIdat,1),nanstd(BIdat,[],1)/sqrt(size(BIdat,1)),...
    'LineColor',burst_color,'ShadeColor',burst_color); hold on % excitation

p2 =errorshade(pltim,nanmean(nonBIdat,1),nanstd(nonBIdat,[],1)/sqrt(size(nonBIdat,1)),...
    'LineColor',noburst_color,'ShadeColor',noburst_color); hold on % excitation

yL = ylim;
cut_win = [-0.01 0.01];
line([cut_win; cut_win],repmat(yL',[1 2]),'Color','k','LineWidth',1,'LineStyle','--')



xlabel('Counts'); ylabel('Time lag (s)');

txt1 = text(pltim(end),yL(2)*0.2,['BI>=' num2str(BIthr) ', n = ' num2str(size(BIdat,1))],'Color',burst_color,'HorizontalAlignment','right')
txt2 = text(pltim(end),yL(2)*0.1,['BI<' num2str(BIthr) ', n = ' num2str(size(nonBIdat,1))],'Color',noburst_color,'HorizontalAlignment','right')


fnm = fullfile(resdir,['ACGavg_BI_nonBI']);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);


end



%--------------------------------------------------------------------------
function bursting_RT(resdir,param,avgBI,plotside)



global rootdir figdir_pd

if ~isfolder(resdir); mkdir(resdir); end;

SessList = getdata2analyse(rootdir, 'rectype','LFP','rectime','intraop', ...
    'patients','allpatients', 'side',plotside, 'condition','stimoff');

if strcmp(param,'RT')
    
    
    
    [RT,Go_RT,FAlarm_RT,Hit,NoStopTrials StopTrials,perf,stopperf] = calc_RT(SessList,'TrialEvents_nosync_stimoff.mat',false);
    RT2 = cellfun(@(x) nanmedian(rmoutliers(x)),RT,'UniformOutput',false);
    
    
    
    
    cellids = findcell;
    brst = get_prop('bursting',cellids);
    
    if ~avgBI
        rtval = nan(length(brst),1);
        brstval = nan(length(brst),1);;
    else
        sort_brst = cell(1,length(RT2));
        rtval = [RT2{:}];
    end
    for ic = 1:length(cellids)
        [patnm, sessnm, chan, ~] = cellid2tags(cellids{ic});
        if strcmp(sessnm(end),'l'); side = 'left'; else; side = 'right'; end;
        sinx = find(ismember({SessList.side},side));
        pinx = find(ismember({SessList.patient},patnm));
        minx = intersect(sinx,pinx);
        if ~isempty(minx)
            if ~avgBI
                rtval(ic) = RT2{minx};
                brstval(ic) = brst(ic);
            else
                if ~isempty(RT2{minx})
                    sort_brst{minx} = [sort_brst{minx} brst(ic)];
                end
                
                
            end
        end
    end
    if avgBI
        brstval = cellfun(@nanmedian,sort_brst);
    end
    
    
else
    
    
    ssrt_dir = fullfile(figdir_pd,'Behav','SSRT');
    load(fullfile(ssrt_dir,'intraop_stimoff','SSRT_results.mat'));
    
    pats = fieldnames(SSRT_results);
    
    
    
    cellids = findcell;
    
    brst = get_prop('bursting',cellids);
    
    if ~avgBI
        rtval = nan(length(brst),1);
        brstval = nan(length(brst),1);
    else
        [RT_sort,sort_brst] = deal(cell(1,size(SessList,1)));
    end
    for ic = 1:length(cellids)
        [patnm, sessnm, chan, ~] = cellid2tags(cellids{ic});
        
        if strcmp(sessnm(end),'l'); side = 'left'; else; side = 'right'; end;
        sinx = find(ismember({SessList.side},side));
        pinx = find(ismember({SessList.patient},patnm));
        minx = intersect(sinx,pinx);
        
        if ~isempty(minx)
            if ~avgBI
                rtval(ic) = SSRT_results.(patnm).(['S' sessnm]).(param);
                brstval(ic) = brst(ic);
            else
                
                
                
                RT_sort{minx} = SSRT_results.(patnm).(['S' sessnm]).(param);
                if ~isempty(RT_sort{minx})
                    sort_brst{minx} = [sort_brst{minx} brst(ic)];
                end
                
                
            end
        end
    end
    if avgBI
        rtval = [RT_sort{:}];
        brstval = cellfun(@nanmedian,sort_brst);
    end
    
end
nn = ~isnan(brstval);


fig = figure;

%  MED OFF
[p, R] = polypredcicall_mod(brstval(nn),rtval(nn),0.95,'robust',0.1);
hold on;
scatter(brstval(nn),rtval(nn),[],'k','filled');
ylabel('Bursting'); xlabel(param);



fnm = fullfile(resdir,['bursting_' param '_medBI' char(string(avgBI)) '_' plotside ]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);


end




%-------------------------------------------------------------------------
function burst_mediansplit(splitthr,respevents,resdir,partition,alphas)

close all
global group_dir
% Bursting prop
prop = 'bursting';
sigma = 0.08;
% Resp cells
load(fullfile(group_dir,'RespCells.mat'));
% splitthr = 'median'; % zero
burst_tags = {'Higher burst index','Lower burst index'};
wn = [-3 3];
saveplot_wn = [-1.2 3];
time = (wn(1):1/1000:wn(2));
alphnr = length(alphas);
%%

if length(respevents)>1
    pdcellids = cell(1,2);
    for ri = 1:length(respevents)
        pdcellids{1} = cat(1,pdcellids{1} ,RespCells.(respevents{ri}).none.Activ);
        pdcellids{2} = cat(1,pdcellids{2},RespCells.(respevents{ri}).none.Inhib);
    end
    event = 'allevents';
else
    event = respevents{1};
    pdcellids{1} = RespCells.(event).none.Activ;
    pdcellids{2} = RespCells.(event).none.Inhib;
    
    switch event
        case 'StimulusOn'
            colors = [0 0 1; 0.3010, 0.7450, 0.9330];
            bwin = [-2.5 -1];
        case 'StopSignal'
            colors = [0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250];
            bwin = [-3 -1.5];
        case 'KeyPress1'
            colors = [0.4940, 0.1840, 0.5560; 0.75, 0, 0.75];
            bwin = [-3 -1.5];
        case 'Feedback'
            colors = [0 1 0; 0.4660, 0.6740, 0.1880 ];
            bwin = [-3 -1.5];
    end
end

resptypes = {'activation','inhibition'};



%% PSTHs
group_ev_dir1 = fullfile(resdir,'PSTH_maps'); if ~isdir(group_ev_dir1); mkdir(group_ev_dir1); end;
group_ev_dir2 = fullfile(resdir,'PSTH_avg'); if ~isdir(group_ev_dir2); mkdir(group_ev_dir2); end;


avgfig = figure(3);
for i = 1:2 % celltypes
    
    
    %% PSTH maps
    fig = figure(2);
    cells = pdcellids{i};
    [bix cellids2] = get_prop(prop,cells);
    
    if strcmp(splitthr,'median')
        med = nanmedian(bix);
        gr{1} = cells(bix>=med); % high-bursting
        gr{2} = cells(bix<med); % low-bursting
        
    elseif isnumeric(splitthr)
        gr{1} = cells(bix>splitthr); % high-bursting
        gr{2} = cells(bix<=splitthr); % low-bursting
    end
    
    if strcmp(partition,'all')
        subplot(2,1,1)
        [psth_R1] = norm_psth_map1(gr{1},resptypes{i},event,'bwin', bwin,...
            'sigma',sigma,'parts',partition,'isfig',true,'cLim',[-10 10]);
        title(burst_tags{1})
        
        subplot(2,1,2)
        [psth_R2]  = norm_psth_map1(gr{2},resptypes{i},event,'bwin', bwin,...
            'sigma',sigma,'parts',partition,'isfig',true,'cLim',[-10 10]);
        title(burst_tags{2})
        
        xL = xlim;
        subplot(2,1,1); xlim(interp1(wn,xL,saveplot_wn));
        subplot(2,1,2); xlim(interp1(wn,xL,saveplot_wn));
        
        set(fig,'Position',get(0,'Screensize'))
        fnma = fullfile(group_ev_dir1,['bursting_' num2str(splitthr) 'split_' event '_' resptypes{i}]);
        saveas(fig,[fnma '.jpg'])
        saveas(fig,[fnma '.fig'])
        saveas(fig,[fnma '.pdf'])
        close(fig)
        
        %% PSTH AVG
        figure(3);
        subplot(2,1,i)
        avg_psth_plot(psth_R1,psth_R2, event, colors)
        ylim([-9 9])
        xlim(saveplot_wn);
        legend(burst_tags,'AutoUpdate','off');
        title(resptypes{i})
        for ai = 1:alphnr
            %          [~, pgroup, ~, ~, ~, ~] = std_stat({psth_R1',psth_R2'},'condstats','off','groupstats','on','mode','eeglab',...
            %             'method','perm','mcorrect','fdr','alpha',0.05,'naccu',1000);
            [~, pgroup, ~,~, ~,~] = std_stat({psth_R1',psth_R2'},'condstats','off','groupstats','on',...
                'mode','fieldtrip','fieldtripmethod','montecarlo','mcorrect','cluster',...
                'fieldtripalpha',alphas(ai),'fieldtripclusterparam',{'clusterstatistic','wcm'});
            
            if ~isempty(pgroup)
                %                 sign = find(pgroup{1});
                sign_p = pgroup{1};
                hold on; yL = ylim;
                %                 scatter(time(sign),ones(1,length(sign))*yL(1),30,'k','filled')
                
                if ai==1
                    draw_signifpatch(time,sign_p,[0.8 0.8 0.8])
                    %             scatter(plottime(sign_p),repmat(yL(1)*0.9,1,sum(sign_p)),[],[0.3 0.3 0.3],'filled')
                elseif ai==2
                    draw_signifpatch(time,sign_p,[0.5 0.5 0.5])
                elseif ai==3
                    draw_signifpatch(time,sign_p,[0.3 0.3 0.3])
                end
                
                
            end
        end
    else
        
        [partlabels, partcolors] = makeColorsLabels(@defineLabelsColors_pd,{[partition(2:end) '=1'],[partition(2:end) '=2']});
        
        for kk = 1:2 % high/ low bursting
            
            fig = figure(kk);
            [psth_R{kk}] = norm_psth_map1(gr{kk},resptypes{i},event,'bwin', bwin,...
                'sigma',sigma,'parts',partition,'isfig',true,'cLim',[-5 5]);
            
            xL = xlim;
            subplot(2,1,1); xlim(interp1(wn,xL,saveplot_wn));
            subplot(2,1,2); xlim(interp1(wn,xL,saveplot_wn));
            suptitle(burst_tags{kk})
            
            psth_part1 = psth_R{kk}{1}; psth_part1(find(isnan(psth_part1(:,1))),:) = [];
            psth_part2 = psth_R{kk}{2}; psth_part2(find(isnan(psth_part2(:,1))),:) = [];
            
            set(fig,'Position',get(0,'Screensize'))
            fnma = fullfile(group_ev_dir1,['2PARTS_bursting_' num2str(splitthr) 'split_' event '_' resptypes{i} '_' burst_tags{kk}]);
            saveas(fig,[fnma '.jpg'])
            saveas(fig,[fnma '.fig'])
            saveas(fig,[fnma '.emf'])
            close(fig)
            
            
            
            figure(3)
            
            subplot(2,1,kk)
            avg_psth_plot(psth_part1,psth_part2, event, cat(1,partcolors{:}))
            ylim([-4 4])
            legend(partlabels,'AutoUpdate','off');
            title(burst_tags{kk})
            xlim(saveplot_wn);
            
            for ai = 1:alphnr
                %             [~, pgroup, ~,~, ~,~] = std_stat({psth_part1',psth_part2'},'condstats','off','groupstats','on','mode','eeglab',...
                %                 'method','perm','mcorrect','fdr','alpha',0.05,'naccu',1000);
                [pcond, ~, ~,~, ~,~] = std_stat({psth_part1';psth_part2'},'condstats','on','groupstats','off',...
                    'mode','fieldtrip','fieldtripmethod','montecarlo','mcorrect','cluster',...
                    'fieldtripalpha',alphas(ai),'fieldtripclusterparam',{'clusterstatistic','wcm'});
                %
                %             [pcond, ~, ~, statscond, ~, ~] = std_stat({psth_R{kk}{1}';psth_R{kk}{2}'},'condstats','on','groupstats','off','mode','eeglab',...
                %                 'method','perm','mcorrect','fdr','alpha',0.05,'naccu',50);
                if ~isempty(pcond)
                    %                 sign = find(pcond{1});
                    hold on;
                    %                     yL = ylim;
                    %                 scatter(time(sign),ones(1,length(sign))*yL(1),50,'k','filled')
                    sign_p = pcond{1};
                    if ai==1
                        draw_signifpatch(time,sign_p,[0.8 0.8 0.8])
                        %             scatter(plottime(sign_p),repmat(yL(1)*0.9,1,sum(sign_p)),[],[0.3 0.3 0.3],'filled')
                    elseif ai==2
                        draw_signifpatch(time,sign_p,[0.5 0.5 0.5])
                    elseif ai==3
                        draw_signifpatch(time,sign_p,[0.3 0.3 0.3])
                    end
                    
                end
            end
        end
        
        set(0, 'DefaultFigureRenderer', 'painters');
        suptitle(resptypes{i})
        %         set(figure(3),'Position',get(0,'Screensize'))
        fnmb = fullfile(group_ev_dir2,['2PARTS_bursting_' num2str(splitthr) 'split_' event '_' resptypes{i}]);
        saveas(figure(3),[fnmb '.jpg'])
        saveas(figure(3),[fnmb '.fig'])
        saveas(figure(3),[fnmb '.pdf'])
        close(figure(3))
        
    end
end

if strcmp(partition,'all')
    set(0, 'DefaultFigureRenderer', 'painters');
    %     set(avgfig,'Position',get(0,'Screensize'))
    fnmb = fullfile(group_ev_dir2,['bursting_' num2str(splitthr) 'split_' event]);
    saveas(figure(3),[fnmb '.jpg'])
    saveas(figure(3),[fnmb '.fig'])
    saveas(figure(3),[fnmb '.pdf'])
    close(figure(3))
end



%% PSTH separate line plots


% group_ev_dir3 = fullfile(resdir,'PSTH_lineplots'); if ~isdir(group_ev_dir3); mkdir(group_ev_dir3); end;
% for i = 1:2 % celltypes
%
%
%     %% PSTH maps
%     cells = pdcellids{i};
%     [bix cellids2] = get_prop(prop,cells);
%
%     if strcmp(splitthr,'median')
%         med = nanmedian(bix);
%         gr{1} = cells(bix>med); % high-bursting
%         gr{2} = cells(bix<=med); % low-bursting
%
%     elseif isnumeric(splitthr)
%         gr{1} = cells(bix>splitthr); % high-bursting
%         gr{2} = cells(bix<=splitthr); % low-bursting
%     end
%
%
%     [psth_R1] = norm_psth_map1(gr{1},resptypes{i},event,'bwin', bwin,...
%         'sigma',sigma,'parts','all','isfig',false);
%     figure; plot(psth_R1')
%     legend(gr{1})
%
%
%     [psth_R2]  = norm_psth_map1(gr{2},resptypes{i},event,'bwin', bwin,...
%         'sigma',sigma,'parts','all','isfig',false);
%     figure; plot(psth_R2')
%     legend(gr{2})
%
%     [mi1, milat1] = min(psth_R1(:,3000:5000),[],2);
%     [mi2, milat2] = min(psth_R2(:,3000:5000),[],2);
%     fig = figure;
%     boxplot([mi1 mi2]); hold on;
%     scatter(ones(1,length(mi1)),milat1,'filled')
%     scatter(ones(1,length(mi2))*2,milat2,'filled');
%     [h,p] = ttest(mi1,mi2);
%     if p<0.05; colp = 'r'; else; colp = 'k'; end;
%
%     annotation('textbox',[0 1 1 0],'string',['p = ' num2str(p)],...
%         'Color',colp,'LineStyle','none');
%
%     fnmb = fullfile(group_ev_dir3,['bursting_' num2str(splitthr) 'split_' event '_' resptypes{i} '_boxpl']);
%     saveas(fig,[fnmb '.jpg'])
%     saveas(fig,[fnmb '.fig'])
%     saveas(fig,[fnmb '.pdf'])
%     close(fig)
% end

end
