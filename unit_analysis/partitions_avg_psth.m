function partitions_avg_psth(pdcells,twinds,alignevent,partition,propname_resp,condition,isplot,varargin)
%PARTITIONS_AVG_PSTH
%   PARTITIONS_AVG_PSTH(pdcells,twinds,alignevent,partition,propname_resp,condition,isplot)
%         Selects units with significant predictive behaviour (see RESPSORT_PARTITIONS_PD) around ALIGNEVENT.
%         Groups units based on their prediction type.
%         Creates average PSTHs, PSTH maps and pie charts of predictive cells with similar properties.
%         Cells on PSTH maps are sorted in descending/ascending order (depending on
%         selection criteria, CONDITION).
%         Saves 'PredCells.mat' file into group_dir folder, containing all
%         predictive cellids + grouped also according to conditions.

%     Inputs:
%     PDCELLS       cell array of unit IDs to analyse
%     
%     TWINDS        nx2 vector test windows in sec (ex. [-0.5 0; -1 0; -1.5 0; -2 0])
%     
%     ALIGNEVENT	character array, label of behavioral event
%     
%     PARTITION     character array, label of partition criteria
%     
%     PROPNAME_RESP name of property in CellBase corresponding to ultimate
%                   psth stat testing for responsive units
%     
%     CONDITION     character array or cell array of selection criteria ("condition") (ex 'bursting', {'bursting','deltarhythmic'})
%     
%     ISPLOT        1|0, if 1 saves figures
%
%
% See also: RESPSORT_PARTITIONS_PD

% Johanna Petra Szabó, Balázs Hangya, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu



dbstop if error
global cell_dir group_dir


if nargin<8
    excl_beforeEv =[];
else
    excl_beforeEv = varargin{1};
end

if nargin<9
    save_PredCells =false;
else
    save_PredCells = varargin{2};
end

numCells = length(pdcells);

propnames = listtag('prop');


try
    load(fullfile(group_dir,'PredCells.mat'));
    
catch
    PredCells = struct;
end



TEtags = {[partition(2:end) '=1'],[partition(2:end) '=2']};
if isempty(excl_beforeEv)
    bwin = [-3 -2];
else
    bwin = [-1.1 -0.1];
end
[parttags, colors] = makeColorsLabels(@defineLabelsColors_pd,TEtags);


% get cells responsive to the event

% propnm = propnames(cellfun(@(x) contains(x,propname_resp),propnames));
% % cell2mat(cellfun(@(x) getvalue(x,pdcells),propnm,'UniformOutput',false));
% resp_values0 = getvalue(propnm{1},pdcells);
% 
% cellidx_act = find(resp_values0==1);
% cellidx_inh = find(resp_values0==-1);
% 
% cellids_act = pdcells(cellidx_act);
% cellids_inh = pdcells(cellidx_inh);

[isresp, cids] = get_prop('resp',pdcells, 'Events2Align',{alignevent});

cellids_act = cids(isresp==1);
cellids_inh = cids(isresp==-1);


% get cells predictive to event partitions

pred_codes = [-1 1]; pred_tags = {'F<S','F>S'};

 propnm_pred= propnames(cellfun(@(x) contains(x,alignevent) & contains(x,partition(2:end)) & contains(x,'_stat_'),propnames));
    
    if isempty(propnm_pred)
        propnm_pred = propnames(cellfun(@(x) contains(x,partition(2:end)) & contains(x,'_stat_') & ~contains(x,'StimulusOn'),propnames));
        disp('StopSignal partitioned')
    end
    
    
    pred_values = cell2mat(cellfun(@(x) getvalue(x,pdcells), propnm_pred,'UniformOutput',false));
    
    
    pred_cells_inx = find(any(pred_values,2));
    
    all_pred_values = pred_values(pred_cells_inx,:); % predictive values (-1/0/1) derived from statistical tests with different test windows
    opmat = mat2cell(all_pred_values,ones(size(all_pred_values,1),1),size(all_pred_values,2));
    
    values = cellfun(@(x) x(find(x==1|x==-1,1)), opmat); % values of first test for each cell
    
    inx1 = values==pred_codes(1);
    inx2 =  values==pred_codes(2);
    
if save_PredCells
   
    pred_cellidx1 = pred_cells_inx(inx1);
    pred_cellidx2 = pred_cells_inx(inx2);
else
    pred_cellids1 = PredCells.([alignevent '_StopPartition']).none.F_smaller_than_S;
    pred_cellids2 = PredCells.([alignevent '_StopPartition']).none.F_bigger_than_S;
    
    pred_cellidx1 = find(ismember(pdcells,pred_cellids1));
    pred_cellidx2 = find(ismember(pdcells,pred_cellids2));
    
end

% plot cells grouped by conditions
cnr = length(condition);
for coi = 1:cnr
    
    cond = condition{coi};
    [pred_cellids11, pred_cellids22, bindex1, bindex2, ~] = find_cellidx(cond,cell_dir,pred_cellidx1,pred_cellidx2,pdcells); % find cells matching condition
    
    
    PredCells.([alignevent '_' partition(2:end)]).(cond).F_smaller_than_S = pred_cellids11;
    PredCells.([alignevent '_' partition(2:end)]).(cond).F_bigger_than_S = pred_cellids22;
    
    try
        PredCells.([alignevent '_' partition(2:end)]).(cond) = ...
            rmfield( PredCells.([alignevent '_' partition(2:end)]).(cond),'Activ');
        PredCells.([alignevent '_' partition(2:end)]).(cond) = ...
            rmfield( PredCells.([alignevent '_' partition(2:end)]).(cond),'Inhib');
    end
    
    if isempty(pred_cellids11) && isempty(pred_cellids22)
        continue
    end
    
    
    if isplot
        % Plot PSTH maps
        cellnrs = [length(pred_cellids11) length(pred_cellids22)];
        all_cellids = cat(2,pred_cellids11,pred_cellids22);
        
        Ha = figure;
        [R] = norm_psth_map1(all_cellids,'grouped',alignevent,...
            'baseline','indiv','basl_psth',{},'bwin',bwin, 'parts',partition,'cLim',[-6 6],...
            'bindex',cat(1,bindex1,bindex2),'isfig',true,'grouptags',pred_tags,'group_limit',cellnrs(1),...
            'excl_beforeEv',excl_beforeEv);
        
        
        % Mark significantly different test window
        allcellnr = sum(cellnrs);
        all_preds_orig = cat(2,pdcells(pred_cellidx1),pdcells(pred_cellidx2));
        condinx= ismember(all_preds_orig,all_cellids);
        
        all_pred_values_sorted = cat(1,all_pred_values(inx1,:),all_pred_values(inx2,:));
        all_pred_values_sorted = all_pred_values_sorted(condinx,:);
        
        all_pred_values_sorted(all_pred_values_sorted~=0) = 1; % =1 if sign. diff (independently from direction)
        
        statint = twinds(:,1)*1000; % intervals of test windows in ms
        
        xL = Ha.Children(2).XLim;
        Xdotinx =  interp1([-3000 3000],xL,statint);
        Xdotinx_rep = repmat(Xdotinx',allcellnr,1);
        
            xmat = all_pred_values_sorted .* Xdotinx_rep;
            ymat = all_pred_values_sorted .* repmat((1:allcellnr)',1,size(twinds,1));
   
        hold on; scatter(xmat(:), ymat(:),30,'red','filled')
        
        
        % Save PSTH maps + matrix
        fdir2 = fullfile(group_dir,alignevent,'pred_PSTH_maps'); fastif(~isdir(fdir2),mkdir(fdir2),0);
        
        fnm3 = fullfile(fdir2,[alignevent '_average.mat']);   % full PSTH matrix
        fnm = fullfile(fdir2,[alignevent '_' partition(2:end) '_' cond]);
        
        set(0, 'DefaultFigureRenderer', 'painters');
        
        if ~isempty(excl_beforeEv)
            fnm3 = [fnm3 '_ExcBef_' excl_beforeEv];
        end
        save(fnm3,'R');
        
        
        if ~isempty(excl_beforeEv)
            fnm = [fnm '_ExcBef_' excl_beforeEv];
        end
        saveas(Ha,[fnm '.jpg']);
        saveas(Ha,[fnm '.fig']);
%         saveas(Ha,[fnm '.pdf']);
        close(Ha)
        
        
        
        %% Plot AVG PSTH
        
        in1 = 1:cellnrs(1);
        in2 = cellnrs(1)+1:allcellnr;
        
        H = figure;
        set(H,'Visible','off');
        subplot(3,4,[1 2])
        avg_psth_plot(R{1}(in1,:),R{2}(in1,:), alignevent, [colors{1}; colors{2}])
        
        title(pred_tags{1})
        
        subplot(3,4,[5 6])
        avg_psth_plot(R{1}(in2,:),R{2}(in2,:), alignevent, [colors{1}; colors{2}])
        
        title(pred_tags{2})
        
        
        leg = legend(parttags,'Location','best');
        set(leg,'AutoUpdate','off');
        
        % Pies
        subplot(3,4,3)
        
        a = intersect(cellids_act,pred_cellids11);
        cellnr_pie([length(a) length(cellids_act)-length(a)],[],{'predict','no predict'});
        title('Activated Stop-resp', 'Units', 'normalized', 'Position', [0.5, 1.8, 0])
        
        subplot(3,4,4)
        
        i = intersect(cellids_inh,pred_cellids11);
        cellnr_pie([length(i) length(cellids_inh)-length(i)],[],{'predict','no predict'});
        title('Inhibited Stop-resp', 'Units', 'normalized', 'Position', [0.5, 1.8, 0])
        
        
        
        
        
        
        subplot(3,4,7)
        
        a = intersect(cellids_act,pred_cellids22);
        cellnr_pie([length(a) length(cellids_act)-length(a)],[],{'predict','no predict'});
        title('Activated Stop-resp', 'Units', 'normalized', 'Position', [0.5, 1.8, 0])
        
        subplot(3,4,8)
        
        i = intersect(cellids_inh,pred_cellids22);
        cellnr_pie([length(i) length(cellids_inh)-length(i)],[],{'predict','no predict'});
        title('Inhibited Stop-resp', 'Units', 'normalized', 'Position', [0.5, 1.8, 0])
        
        
        
        
        subplot(3,4,10)
        
        cellnr_pie([cellnrs(1) cellnrs(2) numCells-(cellnrs(1)+cellnrs(2))],[],{pred_tags{1},pred_tags{2},'no predict'});
        title('All cells', 'Units', 'normalized', 'Position', [0.5, 1.8, 0])
        
        
        set(H,'Position',get(0,'Screensize'));
        
        
        % Save
        fdir1 = fullfile(group_dir,alignevent,'pred_PSTH_avg'); fastif(~isdir(fdir1),mkdir(fdir1),0);
        fnm1 = fullfile(fdir1,[alignevent '_' partition(2:end) '_' cond '_average']);
        
        if ~isempty(excl_beforeEv)
            fnm1 = [fnm1 '_ExcBef_' excl_beforeEv];
        end
        
        set(0, 'DefaultFigureRenderer', 'painters');
        
        saveas(H,[fnm1 '.jpg']);
        saveas(H,[fnm1 '.pdf']);
        saveas(H,[fnm1 '.fig']);
        
        close(H)
    end
    
end

if save_PredCells
save(fullfile(group_dir,'PredCells.mat'),'PredCells');
end




