function PC_respcells_stacked_compare(pcevent,plot_win,fr_name,downsamp,PC_win,alpha,rectype,chanmean,subregion,PCdir)
% PC_RESPCELLS_STACKED_COMPARE
%   PC_RESPCELLS_STACKED_COMPARE(EventTypes,pcevent,plot_win,fr_name,downsamp,PC_win,alpha,rectype,chanmean,subregion,PCdir)
%       Plots stacked bar plots with count of significantly coupled units
%       within unit groups (responsive/ predicive units). 
% Input parameters:
%     PCEVENT           char. array, label of event around which
%                       spike-phase coupling is calculated
%     PLOT_WIN          1x2 vector, time window relative to event timestamp in
%                       sec, to use for SPC calculation
%     FR_NAME           char. array, label of frequency band of interest
% 
%     DOWNSAMP          spike or trial nr. downsampled
%                    'no' | 'spike' 
%     PC_WIN            nr of smaller time windows (time window of plot_win will
%                       be divided to PC_WIN nr of smaller windows to use for
%                       SPC calculation) 
%     ALPHA             numeric, alpha level, above this threshold units
%                       are considered sign. coupled)
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


global group_dir

EventTypes = {'StimulusOn','StopSignal','KeyPress1','Feedback'};


% Resp cells
load(fullfile(group_dir,'RespCells.mat'));
load(fullfile(group_dir,'PredCells.mat'));

lab = {};
brst = {};
pdcellids = cell(length(EventTypes),2);
for ei = 1:length(EventTypes)
    respevent = EventTypes{ei};
    
    
    pdcellids{ei,1} = RespCells.(respevent).none.Activ;
    lab = cat(2,lab,{[respevent ' A']});
    pdcellids{ei,2} =RespCells.(respevent).none.Inhib;
    lab = cat(2,lab,{[respevent ' I']});
    if strcmp(respevent,'StopSignal')
        pdcellids{length(EventTypes)+1,1} = PredCells.([respevent '_StopPartition']).none.FailedStop;
        lab = cat(2,lab,{[respevent ' F<S']});
        pdcellids{length(EventTypes)+1,2} = PredCells.([respevent '_StopPartition']).none.SuccesfulStop;
        lab = cat(2,lab,{[respevent ' F>S']});
    end
    
    
end


% PC
if contains(rectype,'LFP')&& chanmean
    chtit = ['chanmean_' subregion 'STN'];
elseif contains(rectype,'LFP')&& ~chanmean
    chtit = 'by-channel';
elseif strcmp(rectype,'EEG')
    chtit = 'F4';
end



resdir = fullfile(PCdir,chtit,pcevent,pcevent,[fr_name '_' num2str(plot_win)]);
load(fullfile(resdir,['PC_results_ds' downsamp '_' num2str(PC_win) 'win.mat']));

allcellids = fieldnames(PC_results.Hilb_PC);
if PC_win==1
    ray_sign = structfun(@(x) x.Ray_P<=alpha,PC_results.Hilb_PC);
    rv = structfun(@(x) x.ResVect,PC_results.Hilb_PC);
    ftms = angle(rv); ftms_sign = ftms(ray_sign);
    mrls = abs(rv); mrls_sign = mrls(ray_sign);
elseif C_win==2
    ray_sign1 = structfun(@(x) x.Ray_P(1)<=alpha,PC_results.Hilb_PC);
    ray_sign2 = structfun(@(x) x.Ray_P(2)<=alpha,PC_results.Hilb_PC);
    ray_sign = ray_sign1|ray_sign2;
end

for ic = 1:length(allcellids)
    allcellids{ic}(end-1) = '.';
end


cellinxx = cellfun(@(x) find(ismember(allcellids,x)),pdcellids,'UniformOutput',0);
cellinxx_sign = cellfun(@(x) find(ismember(allcellids(ray_sign),x)),pdcellids,'UniformOutput',0);

respsign = cellfun(@(x) sum(ray_sign(x)),cellinxx)';
respnum = cellfun(@length,pdcellids)';
respnosign = respnum-respsign;

%% Stacked bar with nr of sign. PC cells
actty = {'Cue A','Cue I', 'Stop A', 'Stop I', 'KP A', 'KP I', 'F A', 'F I', 'Stop F<S', 'Stop F>S'};
fig = figure;


barmat = cat(2,respsign(:),respnosign(:))
ba = bar(barmat,'stacked','FaceColor','flat');

ba(1).CData = [0.5 0.5 0.5];
ba(2).CData = [0 0 0];

xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', actty)
legend({'Sign PC','Not sign. PC'})

set(gcf,'Position',get(0,'Screensize'));

fnm = fullfile(resdir,'signPC_respcells_stacked');
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.emf'])
close(fig);

% %% Compare
% resdir2 = fullfile(resdir,'compare_respcells');
% if ~isfolder(resdir2); mkdir(resdir2); end;
% 
% cellinxx2 = reshape(cellinxx_sign',[size(cellinxx_sign,1)*size(cellinxx_sign,2) 1]);
% ctynr = length(cellinxx2);
% 
% % Compare FTM
% ftmBC = cellfun(@(x) ftms_sign(x), cellinxx2,'UniformOutput',0);
% bnnr = 24;
% fig = figure;
% jjj=1;
% for jj = 1:2:ctynr
%     subplot(1,ctynr/2,jjj)
%     
%     [~, ~, ray_p1,~,~, ftmm1] = PC_polar(ftmBC{jj}',bnnr,'ray',true,false,[.7 .8 1]);
%     
%     hold on;
%     [~,~, ray_p2,~,~, ftmm2] = PC_polar(ftmBC{jj+1}',bnnr,'ray',true,false,[0 0.1 0.2])
%     
%     
%     rL = rlim;
%     polarplot([ftmm1 ftmm1],rL,'Color',[.7 .8 1],'LineWidth',1)
%     polarplot([ftmm2 ftmm2],rL,'Color',[0 0.1 0.2],'LineWidth',1)
%     
%     
%     
%     if ray_p1<0.05; c= 'r'; else; c = 'k'; end;
%     text(pi*1.2,rL(2)*2,['A: p = ' num2str(ray_p1)],'Color',c);
%     
%     if ray_p2<0.05; c= 'r'; else; c = 'k'; end;
%     text(pi*1.25,rL(2)*2,['I: p = ' num2str(ray_p2)],'Color',c);
%     
%     if length(ftmBC{jj})>=18 && length(ftmBC{jj+1})>=18
%         [~, w2_p] = b_watsontwo(ftmBC{jj},ftmBC{jj+1});
%         if w2_p(1)<0.05; c= 'r'; else; c = 'k'; end;
%         text(pi*1.3,rL(2)*2,['A vs I: p=[' num2str(w2_p(1)) ' ' num2str(w2_p(2)) ']'],'Color',c);
%     end
%     legend(actty{jj:jj+1},'Location','northoutside')
%     jjj =jjj+1;
% end
% set(gcf,'Position',get(0,'Screensize'));
% fnm = fullfile(resdir2,'signPC_respcells_FTM');
% saveas(fig,[fnm '.jpg'])
% saveas(fig,[fnm '.fig'])
% saveas(fig,[fnm '.emf'])
% close(fig);
% 
% 
% % Compare MRL
% 
% mrlBC = cellfun(@(x) mrls_sign(x), cellinxx2,'UniformOutput',0);
% mL= max(cellfun(@length,mrlBC));
% mrlBC2 = cellfun(@(x) [x; nan(mL-length(x),1)], mrlBC, 'UniformOutput',0);
% mrl_bar = cat(2,mrlBC2{:});
% 
% fig = figure;
% 
% boxplot(mrl_bar,actty);
% setmyplot_balazs(gca);
% set_my_boxplot(gca);
% 
% yL = ylim;
% 
% [kwp,anovatab,stats] = kruskalwallis(mrl_bar,actty,'off');
% if kwp<=.05; col = 'r'; else; col = 'k'; end;
% text(1,yL(1)*.9,['Kruskal-Wallis p=' num2str(kwp)],'Color',col);
% 
% 
% for j = 1:2:ctynr
%     [rsp, ~, ~] = ranksum(mrlBC{j},mrlBC{j+1});
%     
%     if rsp<=.05; col = 'r'; else; col = 'k'; end;
%     text(j,yL(2)*.9,['A vs I (Ranksum) p=' num2str(rsp)],'Color',col);
%     
% end
% set(gcf,'Position',get(0,'Screensize'));
% fnm = fullfile(resdir2,'signPC_respcells_MRL');
% saveas(fig,[fnm '.jpg'])
% saveas(fig,[fnm '.fig'])
% saveas(fig,[fnm '.emf'])
% close(fig);




