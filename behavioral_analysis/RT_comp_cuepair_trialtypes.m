function RT_comp_cuepair_trialtypes(conditions)
% RT_COMP_CUEPAIR_TRIALTYPES Compares RT values in high and low conflict trials
%   RT_comp_cuepair_trialtypes(conditions) draws boxplots of RT values 
%   associated with high and low conflict trials for each task condition. 
%   Values are compared with Mann-Whitney U test.
%   High conflict trials: cue (go signal): 3-1 sequence.
%   Low conflict trials: cue (go signal): 1-2 sequence.
%
% Input parameters:
%       CONDITIONS      Nx2 cell array of task conditions to compare
%               -first column corresponds to recording time, second column to
%               DBS stimulation condition
%               -ex.:{'preop','stimoff';'intraop','stimoff';'postop','stimoff';'postop','stimon'};
%
% Johanna Petra Szabó, 10.2025
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


indivfig = true;
sides = {'left','right'};

global filesdir figdir_pd

resdir = fullfile(figdir_pd,'Behav','CuePairs');
if ~isdir(resdir); mkdir(resdir); end


for coci = 1:length(conditions); % loop over conditions
    rectime = conditions{coci,1};
    condi = conditions{coci,2};
    
    sess2analyse = getdata2analyse(filesdir, 'rectype','BEHAV',...
        'rectime',rectime,'side','bothside', 'condition',condi);
    
    sessnr = length(sess2analyse);
    
    
    TEnm = ['TrialEvents_nosync_' condi '.mat'];
    
    [RTs0,~,~,~,~,~,~] = calc_RT(sess2analyse,TEnm,false);
    [~,rmx] = cellfun(@(x) rmoutliers(x),RTs0,'UniformOutput',0);
    
    RTs = RTs0;
    for s = 1:sessnr
        RTs{s}(rmx{s})= NaN;
    end
    
    for s = 1:sessnr
        
        sessdir = sess2analyse(s).sessfolder;
        RT_pat = RTs{s};
        patnm = sess2analyse(s).patient;
        side = sess2analyse(s).side;
        
        
        RTmeds{s} = comp_RT_onepatient(RT_pat,sessdir,condi, indivfig,rectime,patnm,side, resdir);
    end
    
    
    % Patient averages
    
    grnm = {'RTdecrease','RTincrease'};
    patgroups = clinical_groups(grnm);
    
    
    
    for ss = 1:3
        if ss<3
            sid = sides{ss};
            sx = find( ismember({sess2analyse.side},sid) );
        else
            sid = 'bothside';
            sx = 1:length(sess2analyse);
        end
        fig = figure;
%                 for p = 1:2
%                     px = find( ismember({sess2analyse.patient},patgroups{p}) );
%                     xx = intersect(sx,px);
%                     bpmat = cat(1,RTmeds{xx});
%                     subplot(1,2,p)
%                     bp_fig(bpmat);
%                     title(grnm{p})
%                 end
%                     suptitle([rectime ' ' condi ' ' sides{ss}])
%         
        bpmat = cat(1,RTmeds{sx});
        bp_fig(bpmat); hold on;
        
        cg = {'k','b'};
        for p = 1:2
            px = find( ismember({sess2analyse.patient},patgroups{p}) );
            xx = intersect(sx,px);
            Lmat = cat(1,RTmeds{xx});
            Ln{p} = line(repmat([1;2],[1 length(xx)]),Lmat','Color',cg{p},'LineWidth',2);
        end
        legend([Ln{1}(1) Ln{2}(1)], grnm)
        title([rectime ' ' condi ' ' sid])
        
        
        
        figdir = fullfile(resdir, 'Pat_avgs');
        if ~isdir(figdir); mkdir(figdir); end;
        
        saveas(fig, fullfile(figdir,[rectime '_' condi '_' sid '.jpg']))
        saveas(fig, fullfile(figdir,[rectime '_' condi '_' sid '.fig']))
        saveas(fig, fullfile(figdir,[rectime '_' condi '_' sid '.pdf']))
        close(fig)
    end
end
end


%--------------------------------------------------------------------------
function RTmeds = comp_RT_onepatient(RT_pat,sessdir,condi, indivfig,rectime,patnm,side, resdir)

RTmeds = [];
figdir = fullfile(resdir, 'Indiv_pats',[rectime '_' condi]);
if ~isdir(figdir); mkdir(figdir); end;


TEfnm = ['TrialEvents_nosync_' condi '.mat'];
try
    [ord,~,~,~, rev_skip] = cuepair_trialtypes(sessdir,TEfnm);
catch
    fprintf('no TE %s %s\n', patnm, side)
    return;
end

% if isempty(ord)||isempty(rev_skip)
%     fprintf('no revskip %s %s\n', patnm, side)
%     return;
% end
% easyhard{1} = RT_pat(ord); easyhard{2} = RT_pat(rev_skip);

if isempty(ord)||isempty(rev_skip)
    fprintf('no revskip %s %s\n', patnm, side)
    return;
end
easyhard{1} = RT_pat(ord); easyhard{2} = RT_pat(rev_skip);



RTmeds = cellfun(@nanmedian, easyhard);

if indivfig
    mL = max( cellfun(@length,easyhard) );
    bp = cellfun(@(x) [x nan(1,mL-length(x))] ,easyhard,'UniformOutput' ,false);
    bpmat = cat(1,bp{:})';
    fig = figure;
    bp_fig(bpmat);
    
    tit = [rectime ' ' condi ' ' patnm ' ' side ];
    
    title(tit);
    saveas(fig, fullfile(figdir,[patnm '_' side '.jpg']))
    saveas(fig, fullfile(figdir,[patnm '_' side  '.fig']))
    close(fig)
end

end

%--------------------------------------------------------------------------
function bp_fig(bpmat)


boxplot(bpmat,{'Ord','RevSkip'});
set_my_boxplot(gca); setmyplot_balazs(gca);

[pval,~,stats] = ranksum(bpmat(:,1), bpmat(:,2));
yL = ylim; if pval<0.05; col = 'r'; else; col = 'k'; end;
text(0.5,yL(2)*.9,['p = ' num2str(pval)],'Color',col);

ylabel('RT');
end