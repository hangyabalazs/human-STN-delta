function unit_bursting(EventTypes,partition,parttags)
%UNIT_BURSTING      Analysis of bursting properties of units
%   UNIT_BURSTING(EventTypes) 
%       -generates an auto-correlogram map (with units
%       sorted based on bursting indeces (BI) ), averaged auto-correlograms and
%       a pie chart of bursting and non-bursting units. An empiric threshold of
%       0.35 BI is used. 
%       -correlates BI with UPDRS scores and behavioral parameters (RT,
%       SSDp0.5)
%       -separates units to low- and high bursting groups based on median BI
%       
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



%%
% Auto-correlogram map, sorted based on bursting index of units
ACGmap_sorted(BIthr,plot_win)


% Average auto-correlograms of bursting/ non-bursting units
avgACG_BI_nonBI(BIthr,plot_win,resdir,burst_color,noburst_color)


% Percentage of bursting units (pie chart)

bursting_pies(BIthr,colors,resdir,0)
%%
bursting_pies(BIthr,colors,resdir,1)
%%

% Correlation of bursting index with UPDRS scores
figdir = fullfile(figdir_pd,'Intraop_SP','Bursting','UPDRS_corrs'); if ~isfolder(figdir); mkdir(figdir); end;
updrs_bursting_corr(false, figdir)


% Correlation of bursting index with behavioral parameterss
for k = 1:2
    switch k; case 1; param = 'RT'; case 2; param = 'ssd05'; end;
    
    resdir_behav = fullfile(figdir_pd,'Intraop_SP','Bursting','behav_corrs'); if ~isdir(resdir_behav); mkdir(resdir_behav); end;
    
    bursting_RT(resdir_behav,param,false,'bothside')
end

end



%--------------------------------------------------------------------------
function bursting_pies(thr,colors,resdir,suamua)

pdcells = findcell;
brst = get_prop('bursting',pdcells);

fig = figure;

bc = brst>=thr; nbc = brst<thr; nn = isnan(brst);

if ~suamua
    pie([sum(bc) sum(nbc) sum(nn)]);
    ax = gca;
    ax.Colormap = [colors; 0 0 0];
    legend({['BI>=' num2str(thr)],['BI<' num2str(thr)],'Not enough spikes'})

    fnm = fullfile(resdir,['pie_' num2str(thr)]);
else
    props = get_prop('SUA',pdcells);
    
    bc_sua = sum(bc+props'==2);
    bc_mua = sum(bc+(~props')==2);
    nbc_sua = sum(nbc+props'==2);
    nbc_mua = sum(nbc+(~props')==2);
    nn_sua = sum(nn+props'==2);
    nn_mua = sum(nn+(~props')==2);
    
    pie([bc_sua bc_mua nbc_sua nbc_mua nn_sua nn_mua]);
    hold on;
    legend({['Bursting SUA, n = ' num2str(bc_sua)], ['Bursting MUA, n = ' num2str(bc_mua)],...
        ['Non-bursting SUA, n = ' num2str(nbc_sua)], ['Non-bursting MUA, n = ' num2str(nbc_mua)], ...
       ['No BI calc. SUA, n = ' num2str(nn_sua)], ['No BI calc. MUA, n = ' num2str(nn_mua)]},'Location','best');
    
    fnm = fullfile(resdir,['pie_' num2str(thr) '_suamua']);

end
    
saveas(fig,[ fnm '.jpg' ])
saveas(fig,[ fnm '.fig' ])
close(fig)


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




