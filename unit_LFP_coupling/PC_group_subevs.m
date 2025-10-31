function PC_group_subevs(plot_win,freqs,fr_names,PCdir,EventTypes,SubEventTypes,varargin)
%PC_GROUP_SUBEVS    Compare spike-phase coupling of subevents
%  PC_GROUP_SUBEVS(plot_win,freqs,fr_names,PCdir,EventTypes,SubEventTypes,...)
%  -Compares phase distribution (stat. test: Watson's two sample test of homogeneity)
%   and mean resultant length (MRL, stat. test: Signed rank test) of
%   SUBEVENTTYPES, calculated around EVENTTYPES.
%
% Required inputs:
%     PLOT_WIN          1x2 vector, time window relative to event timestamp in
%                       sec, to use for SPC calculation
%
%     FREQS             Nx2 matrix, boundaries of frequency bands of
%                       interest (N is number of frequency bands to analyse, first column is
%                       the lower limit, second column is the lower limit)
%
%     FR_NAMES          cell array, labels of frequency bands of interest
%                       (nr of cells has to correspond to nr of rows in FREQS)
%
%     PCDIR             path to save results and figures
%
%     EVENTTYPES        1xN cell array of event labels, ex: {'StimulusOn','StopSignal','KeyPress1','Feedback'};
%
%     SUBEVENTTYPES     Nx2 cell array of partition ("subevent") labels, each row
%                       corresponds to an event label, columns to partitions
%                       {'CueStim','StopStim';'FailedStopTrial','SuccesfulStopTrial';'CueResponse','StopResponse';'Correct','Error';};
%
% Optional input (name-value pairs with default values):
%   'chanmean'      if true, channel averaged LFP data is used, if false
%                   LFP data is plotted channel-by-channel
%                   (relevant only for intraop LFP data) (def. value:
%                   false)
%                   true | false (default value: true)
%   'subregion'     character array or cell array, uses channel data
%                   derived from listed STN subregions
%                   (relevant only for intraop LFP data)
%          'all' | 'Motor' | 'Limbic' | 'Associative'  (default value: 'all')
%   'isplot'        true | false, if 1 figures are generated and saved (def.
%                   value: true)
%   'downsamp'      spike or trial nr. downsampled
%                    'no' | 'spike' (default value: 'spike')
%   'dominantfreq'  if true, dominant frequency within predefined frequency
%                   band limits (FREQS) are used for PC calculation, see FIND_DOMINANT_FREQ_BANDS
%                   true | false (default value: true)
%   'PCwin'         nr of smaller time windows (time window of plot_win will
%                   be divided to PC_WIN nr of smaller windows to use for
%                   SPC calculation) (def. value: 1)
%   'group'         cell array of unit groups;
%        'all' | 'signPC' | [EVENT ' resp'] | [EVENT ' pred'], where EVENT
%        is a char. array of an event label
%   'plottype'      char. array, type of polar plot to represent phase
%                   distribution
%       'scatter' | 'hist' (def. value: 'scatter')
%   'partition'     char array of partition label, for color and
%                   contingency label (see defineLabelsColors_pd.m)
%                   ex: #StopPartition | #CuepairPartition (def: '#StopPartition')
%   'parttags'      vector of partition tag (see defineLabelsColors_pd.m)
%
%
% See also: PC_CELL_LEVEL, PC_GROUPS, B_RAO3_MOD, B_WATSONTWO

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

prs = inputParser;
addRequired(prs,'plot_win',@ismatrix);
addRequired(prs,'freqs',@ismatrix);
addRequired(prs,'fr_names',@iscell);
addRequired(prs,'PCdir',@isdir);
addRequired(prs,'EventTypes',@iscell);
addRequired(prs,'SubEventTypes',@iscell);
addParameter(prs,'chanmean',false,@islogical);
addParameter(prs,'subregion','all',@ischar);
addParameter(prs,'isplot',true,@islogical);
addParameter(prs,'downsamp','spike',@ischar);
addParameter(prs,'dominantfreq',true,@islogical);
addParameter(prs,'PC_win',1,@isvector);
addParameter(prs,'group','all',@ischar);
addParameter(prs,'plottype','scatter',@ischar);
addParameter(prs,'partition','#StopPartition',@ischar);
addParameter(prs,'parttags',[1 2],@isvector);
parse(prs,plot_win,freqs,fr_names,PCdir,EventTypes,SubEventTypes,varargin{:})
pr = prs.Results;



[~,fold] = fileparts(PCdir);
frnr = length(fr_names);


if contains(fold,'LFP')&& pr.chanmean
    chtit = ['chanmean_' pr.subregion 'STN'];
elseif contains(fold,'LFP')&& ~pr.chanmean
    chtit = 'by-channel';
elseif contains(fold,'EEG')
    chtit = 'F4';
end





%% subevs of predef clusters
for fri = 1:length(fr_names)
    if pr.dominantfreq && ~contains(pr.fr_names{fri},'dom')
        frnm = ['dom_' pr.fr_names{fri}];
    else
        frnm = [pr.fr_names{fri}];
    end
    
    
    
    
    
    grnm = pr.group;
    for ei = 1:length(pr.EventTypes)
        event = pr.EventTypes{ei};
        
        structdir = fullfile(PCdir,chtit,event,event,[frnm '_' num2str(plot_win)],'Groups_modelselection');
        figdir = fullfile(structdir,grnm);
        % Load population SPC measures of cell groups
        load(fullfile(structdir,['PC_groups_ds_' pr.downsamp '_win' num2str(pr.PC_win) '.mat'])) % resultant vector of mean angles corresponding to cells in the actual group
        
        
        
        %%
        grnm2 = grnm; grnm2(isspace(grnm)) = '_';
        %         bestmod = PC_groups.(grnm2).avg_Bestmodel;
        bestmod = 1;
        
        
        for ww = 1:pr.PC_win
            for cp = 1:bestmod(ww)
                if bestmod(ww)>1
                    clusters =  PC_groups.(grnm2).clusters_1;
                    cids = {clusters{1}.(['CL' num2str(cp)]).cellids};
                elseif bestmod(ww)==1
                    cids = PC_groups.(grnm2).cellids;
                end
                for sei = 1:2
                    
                    evty = pr.SubEventTypes{ei,sei};
                    
                    resdir = fullfile(PCdir,chtit,event,evty, [frnm '_' num2str(plot_win)]);
                    
                    
                    load(fullfile(resdir,['PC_results_dsspike_' num2str(pr.PC_win) 'win.mat'])); % resultant vector + stats for each cell
                    
                    
                    
                    resvect = nan(1,length(cids));
                    for ic = 1:length(cids)
                        try
                            resvect(ic) = PC_results.Hilb_PC.(cids{ic}).ResVect(ww);
                        catch
                            resvect(ic) = nan;
                        end
                    end
                    resvect(isnan(resvect)) = [];
                    RV.(evty).(['CL' num2str(cp)]) = resvect;
                end
            end
            
            
            fg = clust_fig_subevents(bestmod(ww),pr.SubEventTypes,ei,RV,pr.plottype,pr.partition,pr.parttags);
            fignm = [ '_ds_spike_' pr.partition(2:end) '_comps' num2str(bestmod(ww)) '_' pr.plottype];
            
            switch pr.PC_win;
                case 1;
                    suptitle({event,num2str(plot_win),[pr.group ' group']});
                    fignm = [pr.group fignm];
                case 2;
                    suptitle({event,['Win' num2str(ww)], [pr.group ' group']});
                    fignm = [pr.group '_WIN' num2str(ww) fignm];
            end
            
            
            set(fg,'Position',get(0,'Screensize'))
            saveas(fg,fullfile(figdir,[fignm '.jpg']));
            saveas(fg,fullfile(figdir,[fignm  '.fig']));
            saveas(fg,fullfile(figdir,[fignm  '.pdf']));
            close(fg);
        end
    end
end
end




%--------------------------------------------------------------------------
function fg = clust_fig_subevents(bestmod,SubEventTypes,ei,RV,plottype,partition,parttags)


fg = figure;
spnr = 1;


switch bestmod
    case 1;
        %         cols{1} = rgb('light red');
        %         cols{2} = rgb('light green');
        
        [~, cols] = makeColorsLabels(@defineLabelsColors_pd,...
            {[partition(2:end) '=' num2str(parttags(1)) ],[partition(2:end) '=' num2str(parttags(2))]});
        
        
    case 2;
        
        cols{1} = rgb('light red'); cols{1} = [cols{1} ; rgb('dark red')];
        cols{2} = rgb('light green'); cols{2} = [cols{2} ; rgb('dark green')];
        
    case 3;
        
        cols{1} = rgb('light red'); cols{1} = [cols{1} ; rgb('dark red')]; cols{1} = [cols{1} ; rgb('magenta')];
        cols{2} = rgb('light green'); cols{2} = [cols{2} ; rgb('dark green')]; cols{2} = [cols{2} ; rgb('olive')];
end




rownr = bestmod; colnr = 2;
for cp = 1:bestmod
    subplot(rownr,colnr,spnr)
    
    wat_p = nan(1,2);  rao_p = nan(1,2); nn = nan(1,2);
    for sei = 1:2
        evty = SubEventTypes{ei,sei};
        cl_rv = RV.(evty).(['CL' num2str(cp)]);
        
        thet = angle(cl_rv);
        rhos = abs(cl_rv);
        
         %
        %         [mu,kappa,~,wp,~] = b_watson(thet);
        %      s
        
        if length(thet)>=4
            [~,rao_p(sei), ~,~ ,mrl,ftm] = b_rao3_mod(thet);
        else
            [~,~, ~,~ ,mrl,ftm] = b_rao3_mod(thet);
        end
        
        if strcmp(plottype,'scatter')
            polim = polarscatter(thet,rhos,50,cols{sei}(cp,:),'filled'); hold on;
            maxrL = 1.1;
            
            PL(sei) = polarplot([angle(ftm) angle(ftm)],[0 mrl],'Color',cols{sei}(cp,:)/2,'LineWidth',3);
        elseif strcmp(plottype,'hist')
            edges=[(0:20:360)/180*pi];
            polim = polarhistogram(thet,edges,'FaceColor',cols{sei}(cp,:)); hold on;
            maxrL = 15;
            
            PL(sei) = polarplot([angle(ftm) angle(ftm)],[0 maxrL],'Color',cols{sei}(cp,:)/2,'LineWidth',3);
        end
        
        rlim([0 maxrL]);
        
        set(gca,'ThetaAxisUnits','radians');
        
       nn(sei) = length(thet);
    end
    title(['Cluster' num2str(cp) ])
%     legend(PL, SubEventTypes,'Location','north')
    
    %     annotation('textbox', [0, 1/cp, 1, 0], 'string', ...
    %         [SubEventTypes{ei,1}(1) ' : W p = ' num2str(wat_p(1))],'LineStyle','none')
    %
    %     annotation('textbox', [0, 0.9/cp, 1, 0], 'string', ...
    %         [SubEventTypes{ei,2}(1) ' : W p = ' num2str(wat_p(2))],'LineStyle','none')
    
    
    if rao_p(1)<0.05; colo = 'r'; else; colo = 'k'; end;
    annotation('textbox', [0, 0.8/cp, 1, 0], 'string', ...
        [SubEventTypes{ei,1}(1) ' : R p = ' num2str(rao_p(1)) ', n=' num2str(nn(1))],'LineStyle','none','Color',colo)
    
    if rao_p(2)<0.05; colo = 'r'; else; colo = 'k'; end;
    annotation('textbox', [0, 0.7/cp, 1, 0], 'string', ...
        [SubEventTypes{ei,2}(1) ' : R p = ' num2str(rao_p(2)) ', n=' num2str(nn(2))],'LineStyle','none','Color',colo)
    
    
    
    rv1 = angle(RV.(SubEventTypes{ei,1}).(['CL' num2str(cp)]))';
    rv2 = angle(RV.(SubEventTypes{ei,2}).(['CL' num2str(cp)]))';
    try
        [u2, w2_p] = b_watsontwo(rv1,rv2);
    catch
        fprintf('Sample size too small for Watson\n');
        w2_p = [NaN NaN];
    end
    
    if w2_p(2)<0.05; colo = 'r'; else; colo = 'k'; end;
    annotation('textbox', [0, 1/cp, 1, 0], 'string', ...
        [SubEventTypes{ei,1}(1) 'vs' SubEventTypes{ei,2}(1) ' : W2 p = ' num2str(w2_p(2))],'LineStyle','none','Color',colo)
    
    
    
    
    spnr = spnr+1;
    subplot(rownr,colnr,spnr);
    
    mrls1 = abs(RV.(SubEventTypes{ei,1}).(['CL' num2str(cp)]))';
    mrls2 = abs(RV.(SubEventTypes{ei,2}).(['CL' num2str(cp)]))';
    
    boxplot([mrls1 mrls2],SubEventTypes(ei,:));
    
    set_my_boxplot(gca)
    
    [wix_p, ~, ~] = signrank(mrls1,mrls2);
    hold on;
    yL = ylim; xL = xlim;
    if wix_p<0.05; colo = 'r'; else; colo = 'k'; end;
    text(xL(1),yL(end),['WilcX p = ' num2str(wix_p)],'Color',colo);
    spnr = spnr+1;
    setmyplot_balazs(gca);
    
end
end