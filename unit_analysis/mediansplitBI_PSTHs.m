function mediansplitBI_PSTHs(EventTypes,partition,parttags,group,patgr)
%MEDIANSPLITBI_PSTHS Peri-event time histogram (PSTH) of units separated by median bursting index
% MEDIANSPLITBI_PSTHS(EventTypes,partition,parttags) 
%       generates average PSTH plots (with standard error) for low- and
%       high bursting groups including all epochs around a selected event/
%       including epochs partitioned based PARTITION and PARTTAGS (see defineLabelsColors_pd.m).
%       The difference between partitioned trials is tested
%       with permutation test with cluster based correction in both groups.
%
%Input parameters:
%     EVENTTYPES        1xN cell array of event label for responsive units 
%
%     PARTITION         char. array of partition label (see defineLabelsColors_pd.m)
%
%     PARTTAGS          tag of partition (see defineLabelsColors_pd.m)
%
% See also: AUTOCORR_PD, DEFINELABELSCOLORS_PD, STD_STAT, NORM_PSTH_MAP1, ULTIMATE_PSTH
%
% Johanna Petra Szabó, 05.2025
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

global figdir_pd
resdir = fullfile(figdir_pd,'Intraop_SP','Bursting'); if ~isdir(resdir); mkdir(resdir); end;



for ei = 1:length(EventTypes)
    
    respev = EventTypes(ei);
       
    burst_mediansplit('median',respev,resdir,partition,parttags,.05,group); % PSTHs for stop trials around RESPEV, partitioned based on outcome of stopping + statistics
end

end

%-------------------------------------------------------------------------
function burst_mediansplit(splitthr,respevents,resdir,partition,parttags,alphas,group)

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

if ~isempty(group)||contains(group,'all')
    
    if strcmpi(group,'SUA')
        props1 = get_prop('SUA',pdcellids{1} );
        pdcellids{1} = pdcellids{1}(props1);
        
        props2 = get_prop('SUA',pdcellids{2} );
        pdcellids{2} = pdcellids{2}(props2);
        
    elseif strcmpi(group,'MUA')
        
        props1 = get_prop('SUA',pdcellids{1} );
        pdcellids{1} = pdcellids{1}(~props1);
        
        props2 = get_prop('SUA',pdcellids{2} );
        pdcellids{2} = pdcellids{2}(~props2);
        
    else
        [patgroups, groups_nm] = clinical_groups({group},'intraop','stimoff','left');
        
        cellids0 = {};
        pats = patgroups{1};
        for j = 1:length(pats)
            patcells = findcell('rat',pats{j});
            cellids0 = cat(2,cellids0,patcells);
        end
        pdcellids0 = pdcellids;
        pdcellids = cellfun(@(x) intersect(x, cellids0),pdcellids0,'UniformOutput',0);
    end
end
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
            'sigma',sigma,'parts',partition,'isfig',true,'cLim',[-10 10],'parttags',parttags);
        title(burst_tags{1})
        
        subplot(2,1,2)
        [psth_R2]  = norm_psth_map1(gr{2},resptypes{i},event,'bwin', bwin,...
            'sigma',sigma,'parts',partition,'isfig',true,'cLim',[-10 10],'parttags',parttags);
        title(burst_tags{2})
        
        xL = xlim;
        subplot(2,1,1); xlim(interp1(wn,xL,saveplot_wn));
        subplot(2,1,2); xlim(interp1(wn,xL,saveplot_wn));
        
        set(fig,'Position',get(0,'Screensize'))
        fnma = fullfile(group_ev_dir1,['bursting_' num2str(splitthr) 'split_' event '_' resptypes{i}]);
        if strcmp(group,'SUA')||strcmp(group,'MUA')
            fnma = [fnma '_' group];
        end
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
                'fieldtripalpha',alphas(ai),'fieldtripclusterparam',{'clusterstatistic','maxsize'});
            
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
                setmyplot_balazs(gca);
                
            end
        end
    else
        
        [partlabels, partcolors] = makeColorsLabels(@defineLabelsColors_pd,...
            {[partition(2:end) '=' num2str(parttags(1))],[partition(2:end) '=' num2str(parttags(2))]});
        
        for kk = 1:2 % high/ low bursting
            
            fig = figure(kk);
            [psth_R{kk}] = norm_psth_map1(gr{kk},resptypes{i},event,'bwin', bwin,...
                'sigma',sigma,'parts',partition,'isfig',true,'cLim',[-5 5],'parttags',parttags);
            
            xL = xlim;
            subplot(2,1,1); xlim(interp1(wn,xL,saveplot_wn));
            subplot(2,1,2); xlim(interp1(wn,xL,saveplot_wn));
            suptitle(burst_tags{kk})
            
            psth_part1 = psth_R{kk}{1}; psth_part1(find(isnan(psth_part1(:,1))),:) = [];
            psth_part2 = psth_R{kk}{2}; psth_part2(find(isnan(psth_part2(:,1))),:) = [];
            
            set(fig,'Position',get(0,'Screensize'))
            fnma = fullfile(group_ev_dir1,[partition(2:end) '_bursting_' num2str(splitthr) 'split_' event '_' resptypes{i} '_' burst_tags{kk}]);
            if strcmp(group,'SUA')||strcmp(group,'MUA')
                fnma = [fnma '_' group];
            end
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
        fnmb = fullfile(group_ev_dir2,[partition(2:end) '_bursting_' num2str(splitthr) 'split_' event '_' resptypes{i}]);
        if strcmp(group,'SUA')||strcmp(group,'MUA')
            fnmb = [fnmb '_' group];
        end
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
    if strcmp(group,'SUA')||strcmp(group,'MUA')
        fnmb = [fnmb '_' group];
    end
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