function avg_psth_stat(event,partition,resdir,plot_win)
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
% See also: PARTITION_TRIALS, RESPONSESORTER_PD, STD_STAT

% Johanna Petra Szabó, Balázs Hangya, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


resptypes = {'Activ','Inhib'};
isfig = true;

time = -3:1/1000:3;
plottime = plot_win(1):1/1000:plot_win(2);
plotinx = dsearchn(time',plottime');

[mylabels, mycolors] = makeColorsLabels(@defineLabelsColors_pd,{[partition(2:end) '=1'],[partition(2:end) '=2']});

 switch event
        case 'StimulusOn'
            bwin = [-2.5 -1];
        case 'StopSignal'
            bwin = [-3 -1.5];
        case 'KeyPress1'
            bwin = [-3 -1.5];
        case 'Feedback'
            bwin = [1.5 3];
 end
    
 %%
 allcells = findcell;

 close all
 cellids = cell(1,length(resptypes));
figure(1);
for rt = 1:length(resptypes)
    
     [isresp, cids] = get_prop('resp',allcells, 'Events2Align',{event});
     
     if strcmp(resptypes{rt},'Activ')
         cellids{rt} = cids(isresp==1);
     elseif  strcmp(resptypes{rt},'Inhib')
          cellids{rt} = cids(isresp==-1);
     end
    
    figure(2)
    psth_n = [];
    [psth_n] = norm_psth_map1(cellids{rt} ,resptypes{rt},event,...
        'baseline','indiv','basl_psth',{},'bwin',bwin, 'parts',partition,'cLim',[-6 6],...
        'bindex',[],'isfig',isfig,'sigma',0.06);
    
    xL = xlim;
    xtime = xL(1):xL(2);
    xplotL= xtime([plotinx(1) plotinx(end)]);
    subplot(211)
    set(gca,'XLim',xplotL)
        subplot(212)
    set(gca,'XLim',xplotL)

    

    fnm =[partition(2:end) '_psth_MAP_' event '_' resptypes{rt}];
    saveas(figure(2),fullfile(resdir,[fnm '.jpg']))
    saveas(figure(2),fullfile(resdir,[fnm '.fig']))
    saveas(figure(2),fullfile(resdir,[fnm '.emf']))
    close(figure(2));
    
    
        figure(1);
    subplot(2,1,rt)
    [pL1,pL2,txt1,txt2] = avg_psth_plot(psth_n{1}(:,plotinx),psth_n{2}(:,plotinx), event, ...
        [mycolors{1}; mycolors{2}],plot_win);
    legend([pL1(end) pL2(end)],mylabels,'AutoUpdate','off','Location','best');
    
    
    statdat = {psth_n{1}';psth_n{2}'};
    
    alphas = [0.05 0.01 0.001];
    for ai = 1:length(alphas)
        [pcond, ~, ~,~,~,~] = std_stat(statdat,'condstats','on',...
            'mode','fieldtrip','fieldtripmethod','montecarlo','mcorrect','cluster',...
            'fieldtripalpha',alphas(ai),'fieldtripclusterparam',{'clusterstatistic','wcm'});
        sign_p = pcond{1};
        
        hold on
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

fnm =[partition(2:end) '_psth_AVG_' event '_resp'];
saveas(figure(1),fullfile(resdir,[fnm '.jpg']))
saveas(figure(1),fullfile(resdir,[fnm '.fig']))
saveas(figure(1),fullfile(resdir,[fnm '.pdf']))
close(figure(1));

end