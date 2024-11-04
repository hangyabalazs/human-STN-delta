function PC_bursting(EventTypes)
% PC_BURSTING   Spike-phase coupling and bursting index
%    PC_BURSTING(event,colors,resdir,chanmean,PC_win)
%       -correlates spike-phase coupling measures with bursting index of
%       significantly coupled units and of responsive units
%       -separates units to low- and high bursting groups based on median bursting index
%       -compares phase distribution (Watson's two sample test of homogeneity) and
%       coupling strength (Mann Whitney U test) low and high bursting units
%       -performs the comparison on all significantly coupled units and on
%       responsive units separately
% Input parameter:
%   EVENTTYPES      cell array of event labels

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


global figdir_pd

burst_color = rgb('maroon');
noburst_color = rgb('gray');
colors = [burst_color; noburst_color];
resdir = fullfile(figdir_pd,'Intraop_SP','Bursting');


% Correlation between spike-phase coupling and bursting index
for eii = 1:length(EventTypes)
    event = EventTypes{eii}; % coupling calculated around EVENT
    
    PC_burst_correlation(event, 'high_delta', false, 'LFP','signPC', true,1, 'no',[-1.5 1.5]);
    
    PC_burst_correlation(event, 'high_delta', false, 'LFP',[event ' resp'], true,1, 'no',[-1.5 1.5]);
end

% Compare phase distribution and coupling strength of low and high bursting, 
%   significantly delta-coupled units
for eii = 1:length(EventTypes) 
    event = EventTypes{eii}; % coupling calculated around EVENT
   
    bursting_mediansplit_PCs_signPC(event,colors,resdir,false,1)
end



% Compare phase distribution and coupling strength of low and high bursting responsive units
for eii = 1:length(EventTypes)
    event = EventTypes{eii}; % coupling calculated around EVENT and selects units responsive to EVENT
    
    bursting_mediansplit_PCs_resp(event,colors,resdir,false,1)
end


end




%--------------------------------------------------------------------------
function bursting_mediansplit_PCs_signPC(event,colors,resdir,chanmean,PC_win)


close all
global cell_dir 
% Bursting prop
prop = 'bursting';

PCdir = fullfile(cell_dir,'PC',['LFP_phase_coupling']);
if chanmean
    chlab = 'chanmean_allSTN';
else
    chlab = 'by-channel';
end
downsamp = 'no';

PCresdir = fullfile(PCdir,chlab,event,event,'dom_high_delta_-1.5         1.5');
load(fullfile(PCresdir,['PC_results_ds' downsamp '_' num2str(PC_win) 'win.mat']))



cellids0 = fieldnames(PC_results.Hilb_PC);
if PC_win==1
    ray_sign = structfun(@(x) x.Ray_P<=0.05,PC_results.Hilb_PC);
elseif PC_win==2
    ray_sign1 = structfun(@(x) x.Ray_P(1)<=0.05,PC_results.Hilb_PC);
    ray_sign2 = structfun(@(x) x.Ray_P(2)<=0.05,PC_results.Hilb_PC);
    ray_sign = ray_sign1|ray_sign2;
end
cells = cellids0(ray_sign);
for ic = 1:length(cells)
    cells{ic}(end-1) = '.';
end





%% PCs

group_ev_dir1 = fullfile(resdir,'mediansplit_PC'); if ~isdir(group_ev_dir1); mkdir(group_ev_dir1); end;




[bix cellids2] = get_prop(prop,cells);

med = nanmedian(bix);
gr{1} = cells(bix>=med);
gr{2} = cells(bix<med);

resvals = cell(1,2);

for j = 1:2 % 1 - higher bursting, 2 - lower bursting
    cellids = gr{j};
    
    resvals{j} = nan(length(cellids),1);
    for ic = 1:length(cellids)
        cellid = cellids{ic};
        cellidF = regexprep(cellid,'\.','_');
        
        try
            resvals{j}(ic) = PC_results.Hilb_PC.(cellidF).ResVect;
        catch
            fprintf('No PC vals %s\n',cellidF); continue;
        end
    end
    
    phasvals = angle(resvals{j});
    polfig = figure(1);
    pl(j) = polarhistogram(phasvals,12,'EdgeColor','k','FaceColor',colors(j,:)); hold on;
    rL = rlim; rmax = rL(2);
    
    
    [~,rao_p, ~,~ ,mrl,ftm] = b_rao3_mod(phasvals');
    
    polarplot([angle(ftm) angle(ftm)],[0 rmax],'Color',colors(j,:),'LineWidth',2)
    
    if rao_p<0.05; colp ='r'; else; colp = 'k'; end;
    annotation('textbox',[0, 1/j, 1,0],'string',['RAO p = ' num2str(rao_p)],'Color',colp,'LineStyle','none');
    
end
legend(pl,{'High burst','Low burst'},'AutoUpdate','off','Location','south outside')
title(event)

if ~any( cellfun(@length,resvals)<18 )
    [u2, w2_p] = b_watsontwo(angle(resvals{1}), angle(resvals{2}));
    if w2_p<0.05; colp ='r'; else; colp = 'k'; end;
    annotation('textbox',[0, 0.1, 1,0],'string',['W2 p = ' num2str(w2_p)],'Color',colp,'LineStyle','none');
end



fnm1 = ['polar_PCs_' event '_signPC'];
saveas(polfig,fullfile(group_ev_dir1,[fnm1 '.jpg']))
saveas(polfig,fullfile(group_ev_dir1,[fnm1 '.fig']))
saveas(polfig,fullfile(group_ev_dir1,[fnm1 '.pdf']))
close(polfig)



mrlfig = figure(2);

maxlen = max(cellfun(@length,resvals));
mrl1 = abs(resvals{1});
mrl2 = abs(resvals{2});
v1 = cat(1,mrl1,nan(maxlen-length(abs(resvals{1})),1));
v2 = cat(1,mrl2,nan(maxlen-length(abs(resvals{2})),1));
boxplot([v1 v2],{'High burst','Low burst'})
hold on;
scatter(   ones(1,length(mrl1)) + randi( [-1 1] ,1,length(mrl1))/10, mrl1,[],'k','filled'  )
scatter(   ones(1,length(mrl2))*2 + randi( [-1 1] ,1,length(mrl2))/10, mrl2,[],'k','filled'  )

[prs,~,~] = ranksum(abs(resvals{1}), abs(resvals{2}))

if prs<0.05; colp ='r'; else; colp = 'k'; end;
yL = ylim;
text(0.5,yL(end)*0.9,['p=' num2str(prs)],'Color',colp);

title(event)



fnm2 = ['mrls_' event '_signPC'];
saveas(mrlfig,fullfile(group_ev_dir1,[fnm2 '.jpg']))
saveas(mrlfig,fullfile(group_ev_dir1,[fnm2 '.fig']))
saveas(mrlfig,fullfile(group_ev_dir1,[fnm2 '.pdf']))
close(mrlfig)

end



%--------------------------------------------------------------------------
function bursting_mediansplit_PCs_resp(event,colors,resdir,chanmean,PC_win)

%
%
%%
close all
global cell_dir group_dir
% Bursting prop
prop = 'bursting';
% Resp cells
load(fullfile(group_dir,'RespCells.mat'));

PCdir = fullfile(cell_dir,'PC',['LFP_phase_coupling']);
if chanmean
    chlab = 'chanmean_allSTN';
else
    chlab = 'by-channel';
end
downsamp = 'no';

PCresdir = fullfile(PCdir,chlab,event,event,'dom_high_delta_-1.5         1.5');
load(fullfile(PCresdir,['PC_results_ds' downsamp '_' num2str(PC_win) 'win.mat']))

pdcellids{1} = RespCells.(event).none.Activ;
pdcellids{2} = RespCells.(event).none.Inhib;

resptypes = {'activation','inhibition'};




%% PCs

group_ev_dir1 = fullfile(resdir,'mediansplit_PC'); if ~isdir(group_ev_dir1); mkdir(group_ev_dir1); end;

polfig = figure(1); mrlfig = figure(2);
for i = 1:2 % 1 - activation, 2 - inhibition
    
    
    cells = pdcellids{i};
    [bix cellids2] = get_prop(prop,cells);
    
    med = nanmedian(bix);
    gr{1} = cells(bix>med);
    gr{2} = cells(bix<=med);
    
    resvals = cell(1,2);
    
    for j = 1:2 % 1 - higher bursting, 2 - lower bursting
        cellids = gr{j};
        
        resvals{j} = nan(length(cellids),1);
        for ic = 1:length(cellids)
            cellid = cellids{ic};
            cellidF = regexprep(cellid,'\.','_');
            
            try
                resvals{j}(ic) = PC_results.Hilb_PC.(cellidF).ResVect;
            catch
                fprintf('No PC vals %s\n',cellidF); continue;
            end
        end
        
        phasvals = angle(resvals{j});
        figure(1); subplot(1,2,i);
        pl(j) = polarhistogram(phasvals,12,'EdgeColor','k','FaceColor',colors(j,:)); hold on;
        rL = rlim; rmax = rL(2);
        
        
        [~,rao_p, ~,~ ,mrl,ftm] = b_rao3_mod(phasvals');
        
        polarplot([angle(ftm) angle(ftm)],[0 rmax],'Color',colors(j,:),'LineWidth',2)
        
        if rao_p<0.05; colp ='r'; else; colp = 'k'; end;
        annotation('textbox',[1-1/i, 1/j, 1,0],'string',['RAO p = ' num2str(rao_p)],'Color',colp,'LineStyle','none');
        
    end
    legend(pl,{'High burst','Low burst'},'AutoUpdate','off','Location','south outside')
    title(resptypes{i})
    
    if ~any( cellfun(@length,resvals)<18 )
        [u2, w2_p] = b_watsontwo(angle(resvals{1}), angle(resvals{2}));
        if w2_p<0.05; colp ='r'; else; colp = 'k'; end;
        annotation('textbox',[0, 0.1, 1,0],'string',['W2 p = ' num2str(w2_p)],'Color',colp,'LineStyle','none');
    end
    
    
    figure(2); subplot(1,2,i);
    
    maxlen = max(cellfun(@length,resvals));
    mrl1 = abs(resvals{1});
    mrl2 = abs(resvals{2});
    v1 = cat(1,mrl1,nan(maxlen-length(abs(resvals{1})),1));
    v2 = cat(1,mrl2,nan(maxlen-length(abs(resvals{2})),1));
    boxplot([v1 v2],{'High burst','Low burst'})
    hold on;
    scatter(   ones(1,length(mrl1)) + randi( [-1 1] ,1,length(mrl1))/10, mrl1,[],'k','filled'  )
    scatter(   ones(1,length(mrl2))*2 + randi( [-1 1] ,1,length(mrl2))/10, mrl2,[],'k','filled'  )
    
    [prs,~,~] = ranksum(abs(resvals{1}), abs(resvals{2}))
    
    if prs<0.05; colp ='r'; else; colp = 'k'; end;
    yL = ylim;
    text(0.5,yL(end)*0.9,['p=' num2str(prs)],'Color',colp);
    
    title(resptypes{i})
    
    
end

figure(1);suptitle(event)
set(figure(1),'Position',get(0,'Screensize'));

figure(2);suptitle(event)

fnm1 = ['polar_PCs_' event];
saveas(figure(1),fullfile(group_ev_dir1,[fnm1 '.jpg']))
saveas(figure(1),fullfile(group_ev_dir1,[fnm1 '.fig']))
saveas(figure(1),fullfile(group_ev_dir1,[fnm1 '.pdf']))
close(figure(1))

fnm2 = ['mrls_' event];
saveas(figure(2),fullfile(group_ev_dir1,[fnm2 '.jpg']))
saveas(figure(2),fullfile(group_ev_dir1,[fnm2 '.fig']))
saveas(figure(2),fullfile(group_ev_dir1,[fnm2 '.pdf']))
close(figure(2))

end



%--------------------------------------------------------------------------

function PC_burst_correlation(event, fr_name, chanmean, rectype,group, dominantfreq,PC_win, ds, plot_win )

global cell_dir

%%
[PC_results_dir, pdmrl, pdphas, pdcellids,frnm,~,resptypes] = get_phas(event, fr_name, chanmean, rectype,group, dominantfreq,PC_win, ds, plot_win);


%% Get bursting values
type_nr = length(pdcellids);

for ty = 1:type_nr
    cellids = pdcellids{ty};
    mrl = pdmrl{ty};
    phas = pdphas{ty};
    phas = mod(phas,2*pi);
    
    if type_nr ==2
        groupnm = [group '_' resptypes{ty}];
    else
        groupnm = group;
    end
    
    
    [brst, ~] = get_prop('bursting',cellids,'Events2Align',{event});
    
    nn = ~isnan(brst);
    fig = figure;
    
    [p, R] = polypredcicall_mod(mrl(nn),brst(nn),0.95,'robust',0.01); hold on;
    xlabel('MRL'); ylabel('BurstIndex');
    frnmtit = frnm; frnmtit(strfind(frnm,'_')) = '-';
    title({'BurstIndex vs PC', ['Cells sign. coupled to ' frnmtit],['Window: +-' num2str(plot_win) ' sec around ' event]})
    
    figdir =  fullfile(cell_dir,'PC','Burst_PC_correlation',[rectype '_PC'],['dom_' fr_name],event); 
    if ~isfolder(figdir); mkdir(figdir); end;
    fnm = fullfile(figdir,['MRL_'  groupnm '_' frnm '_' event '_ds' ds '_' num2str(PC_win) 'win']);

    
    saveas(fig, [fnm '.fig'])
    saveas(fig, [fnm '.jpg'])
    saveas(fig, [fnm '.pdf'])
    close(fig);
    
    
    
    % Phase
    fig = figure;
    
    
    [Rval, pval] = lincirc_corr2(brst,phas);
    phas_phas = cat(1,phas,phas+2*pi);
    bursburs= repmat(brst,[2,1]);
    
    scatter(phas_phas,bursburs,[],'k','filled');
    tcks = [0 pi/2 pi 3*pi/2 2*pi 2*pi+pi/2 3*pi 2*pi+3*pi/2 4*pi];
    xticks(tcks); xticklabels(arrayfun(@num2str,tcks,'UniformOutput',0));
    
    xlim([0 4*pi])
    
    if pval<0.05; col = 'r'; else; col = 'k'; end;
    xL = xlim; yL = ylim;
    text(xL(1),yL(2)*0.9,['p=' num2str(pval) ', R=' num2str(Rval)]);
    
    
    xlabel('Resultant vector'); ylabel('BurstIndex');
    frnmtit = frnm; frnmtit(strfind(frnm,'_')) = '-';
    title({'BurstIndex vs RV', ['Cells sign. coupled to ' frnmtit],['Window: +-2 sec around ' event]})
    
    figdir =  fullfile(cell_dir,'PC','Burst_PC_correlation',[rectype '_PC'],['dom_' fr_name],event); 
    
    if ~isfolder(figdir); mkdir(figdir); end;
    fnm = fullfile(figdir,['RVphase_' groupnm '_' frnm '_' event '_ds' ds '_' num2str(PC_win) 'win']);
    
    saveas(fig, [fnm '.fig'])
    saveas(fig, [fnm '.jpg'])
    saveas(fig, [fnm '.pdf'])
    close(fig);
end
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