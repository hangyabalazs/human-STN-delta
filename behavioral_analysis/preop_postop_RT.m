function preop_postop_RT(grnm,conditions)
%PREOP_POSTOP_RT    RT change relative to preop condition
%   PREOP_POSTOP_RT(grnm,conditions) calculates RT difference between preop
%   and all other task conditions for the specified patient groups. 
%
%   Required inputs:
%       GRNM    cell array of patient group labels
%          {'tremor-dominant','akinetic-rigid','mixed'}, clinical groups
%          {'RTdecrease','RTincrease'}, patient group based on preop-postop RT change
% 
%  CONDITIONS      1x2 cell array of task conditions to compare
%               -first column corresponds to recording time, second column to
%               DBS stimulation condition
%               -ex.:{'preop','stimoff';'intraop','stimoff';'postop','stimoff';'postop','stimon'};

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


global filesdir figdir_pd

resdir = fullfile(figdir_pd,'Behav','preop_postop_RT');
if ~isfolder(resdir); mkdir(resdir); end;

pat_grs = clinical_groups(grnm);
sides = {'left','right'};
grnr = length(pat_grs);
condnr = size(conditions,1);

for si = 1:2
    sid = sides{si};
    
    SessList = {}; RT = {}; medRTs = {};
    for k = 1:grnr
        
        for j = 1:condnr % preop, postop off and on
            TEnm = ['TrialEvents_nosync_' conditions{j,2} '.mat'];
            
            SessList{k,j} = getdata2analyse(filesdir, 'rectype','BEHAV',...
                'rectime',conditions{j,1},'patients', pat_grs{k}, 'side',sid, 'condition',conditions{j,2});
            
            
            [RT{k,j},~,~,~,~,~,~,~] = calc_RT(SessList{k,j},TEnm,false);
            medRTs{k,j} = cellfun(@nanmedian,RT{k,j});
        end
        % if patient nrs are not identical
        if length(unique( cellfun(@length,SessList(k,:)) ))~=1 
            ps1 = {SessList{k,1}.patient};
            for mm = 2:condnr
                ps2 = {SessList{k,mm}.patient};
                ps1 = intersect(ps1,ps2);
            end
            pinx = {};
            for mm = 1:condnr
                ps2 = {SessList{k,mm}.patient};
                [~,~,pinx{mm}] = intersect(ps1,ps2);
            end
            medRTs(k,:) = arrayfun(@(x) medRTs{k,x}(pinx{x}) ,1:condnr,'UniformOutput',0);
        end
    end
    
    
    diffs = {};
    savtit = [sid '_' grnm{:}];
    fig = figure;
    for m = 2:condnr % preop - all other conditions
        
        for k = 1:grnr
            diffs{k} = medRTs{k,m} - medRTs{k,1};
        end
        
        mL = max(cellfun(@length,diffs));
        diffs_nn = cellfun(@(x) [x nan(1,mL-length(x))]'  ,diffs,'UniformOutput',0);
        
        
        subplot(1,condnr-1,m-1)
        
        bp = cell2mat(diffs_nn);
        
        boxplot(bp,grnm);
        set_my_boxplot(gca)
        hold on;
        for k = 1:grnr
            scatter(ones(mL,1)*k + randi([-10 10],mL,1)*0.01 , bp(:,k) ,[],'k','filled' );
        end
        xL = xlim;
        line(xL,[0 0],'Color','k','LineStyle','--');
        title({[conditions{m,1} ' ' conditions{m,2} ' - preop'],[sid ' side']})
        ylabel('median RT (s)');
        
        savtit = [savtit  '_' conditions{m,1} '-' conditions{m,2}];
        setmyplot_balazs(gca)
    end
    
    saveas(fig,fullfile(resdir,[savtit '.jpg']))
    saveas(fig,fullfile(resdir,[savtit '.fig']))
    saveas(fig,fullfile(resdir,[savtit '.pdf']))
    close(fig);
end
end
