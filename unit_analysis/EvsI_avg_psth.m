function EvsI_avg_psth(pdcells,propname,alignevent,condition,isplot)
%EVSI_AVG_PSTH  Average PSTHs, PSTH maps and pie charts
%   EVSI_AVG_PSTH(pdcells,propname,alignevent,condition,isplot)
%     Selects units with significant acitvation/ inhibition (see RESPONSESORTER_PD) to ALIGNEVENT.
%     Groups units based on their response type.
%     Creates average PSTHs, PSTH maps and pie charts of units in the same group. 
%     Units on PSTH maps are sorted in descending/ascending order (depending on
%     selection criteria (CONDITION)).
%     Uses results of statistical testing of responsiveness (see
%     Saves 'RespCells.mat' file into group_dir folder, containing all activated/inhibited responsive
%     cellids + grouped also according to conditions.
% 
%     Inputs:
%     PDCELLS       cell array of unit IDs to analyse
%     
%     PROPNAME      name of property in CellBase corresponding to ultimate psth stat testing
%     
%     ALIGNEVENT	character array, label of behavioral event
%     
%     CONDITION     character array or cell array of selection criteria ("condition") (ex 'bursting', {'bursting','deltarhythmic'})
%     
%     ISPLOT        1|0, if 1 saves figures
%
% See also: RESPONSESORTER_PD

% Johanna Petra Szabó, Balázs Hangya, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


dbstop if error

global cell_dir group_dir
numCells = length(pdcells);

propnames = listtag('prop');


group_ev_dir = fullfile(group_dir,alignevent); fastif(~isdir(group_ev_dir),mkdir(group_ev_dir),0);

try
    load(fullfile(group_dir,'RespCells.mat'));
catch
    RespCells = struct;
end

if iscell(condition)
    cnr = length(condition);
else
    condition = {condition};
    cnr = 1;
end

if ~iscell(propname)
    propname = {propname};
end

switch alignevent
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
        bwin = [1.5 3];
end

% Get results of responsiveness from CellBase
propnm = propnames(cellfun(@(x) (contains(x,alignevent) & contains(x,propname{1})),propnames));
resp_values = cell2mat(cellfun(@(x) getvalue(x,pdcells),propnm,'UniformOutput',false));

cellidx1 = find(resp_values==1); % activated
cellidx2 = find(resp_values==-1); % inhibited cells



% % PSTH of all responsive cells for baseline normalization
% resp_cellids_all = pdcells([cellidx1; cellidx2]);
% wn = [-3 3]; % time window
% dt = 0.001;
% basl_psth = cell(length(resp_cellids_all),1);
% for iC = 1:length(resp_cellids_all)
%     [~, basl_psth{iC},~, ~, ~,~] = ultimate_psth(resp_cellids_all{iC},'trial', alignevent,wn,...
%         'dt',dt,'display',false,'sigma',0.04,'parts','all','isadaptive',0,...
%         'maxtrialno',Inf,'baselinewin',bwin,'testwin',[0 1],'relative_threshold',0.1);
% end



H = figure;
% set(H,'Visible','off');
for coi = 1:cnr
    
    % Select responsive cells matching condition
    cond = condition{coi};
    [cellids1, cellids2, bindex1, bindex2, cond_cellidx] = find_cellidx(cond,cell_dir,cellidx1,cellidx2,pdcells);
    
    cellnrAll = length(cond_cellidx);
    
    RespCells.(alignevent).(cond).Activ = cellids1;
    RespCells.(alignevent).(cond).Inhib = cellids2;
    
    
    
    
    % Plot PSTH maps
    
    if isplot
        resptypes = {'activation','inhibition'};
        for respk = 1:2
            
            
            switch respk;
                case 1
                    cellids_k = cellids1; bindex_k = bindex1; % activ
                case 2
                    cellids_k = cellids2; bindex_k = bindex2; % inhib
            end
            
            
            
            % Normalize PSTHs + draw PSTH maps
            Ha = figure;
            set(Ha,'Visible','off');
            [psth_k] = norm_psth_map1(cellids_k,resptypes{respk},alignevent,...
                'baseline','indiv','basl_psth',{},'bwin',bwin, 'parts','all','cLim',[-6 6],...
                'bindex',bindex_k,'isfig',true);
            
            
            
            
            
            
            % Save PSTH maps
            
            group_ev_dir1 = fullfile(group_ev_dir,'PSTH_maps'); fastif(~isdir(group_ev_dir1),mkdir(group_ev_dir1),0);
            
            switch respk;
                case 1
                    psth_R1 = psth_k;
                    fnm = fullfile(group_ev_dir1,[cond '_' alignevent '_active']);   % sorted normalized PSTHs
                case 2
                    psth_R2 = psth_k;
                    fnm = fullfile(group_ev_dir1,[cond '_' alignevent '_inhib']);   % sorted normalized PSTHs
            end
            
            saveas(Ha,[fnm '.jpg']);
%             saveas(Ha,[fnm '.pdf']);
            saveas(Ha,[fnm '.fig']);
            close(Ha)
        
        end
        
        % Save normalized PSTHs
        fnm3 = fullfile(group_ev_dir1,[cond '_' alignevent '_normPSTHs.mat']);   % full PSTH matrix
        R = [{psth_R1},{psth_R2}];
        save(fnm3,'R');
        
        
        
        %         if length(propname)>1
        %             statint = [500 1000 1500];
        %
        %             [~, all_rvt] = multitest_resp(propname{1,:},propname,pdcells);
        %
        %
        %             xL = Ha.XLim;
        %             Xdotinx =  interp1([-3000 3000],xL,statint);
        %
        %             cf = cellfun(@(x) find(x),mat2cell(all_pred_values_s,...
        %                 size(all_pred_values_s,1),ones(1,size(all_pred_values_s,2))),'UniformOutput',0);
        %
        %             Xs = [];
        %             for twnr = 1:length(Xdotinx)
        %                 Xs = [Xs repmat(Xdotinx(twnr),[1,length(cf{twnr})])];
        %             end
        %             Ys = cat(1,cf{:})';
        %
        %             hold on; scatter(Xs, Ys,8,'red','filled')
        %         end
        
        
        
        
        % Plot AVG PSTH
        H;
        set(H,'Visible','off');
        subplot(cnr+1,2,(coi-1)*2+1)
        avg_psth_plot(psth_R1,psth_R2, alignevent, colors)
        title(cond)
        
        
        % Plot PIE
        subplot(cnr+1,2,(coi-1)*2+2)
        cellnr1 = length(cellids1); cellnr2 = length(cellids2);
        cellnr_pie(cellnr1,cellnr2,cellnrAll-(cellnr1+ cellnr2),colors,{'activation','inhibition','no change'});
        
        
        
    end
    
end
save(fullfile(group_dir,'RespCells.mat'),'RespCells');




% Plot PIE- all cells (ex. rhythmic vs nonrhythmic)
if isplot
    if ~strcmp(cond, 'none')
        
        subplot(cnr+1,2,2*(cnr+1))
        
        cellidx_all = 1:numCells;
        [cellids1 cellids2 bindex1 bindex2 cond_cellidx1] = find_cellidx(condition{1},cell_dir,cellidx_all,cellidx_all,pdcells);
        
        if length(condition)>1
            [cellids1 cellids2 bindex1 bindex2 cond_cellidx2] = find_cellidx(condition{2},cell_dir,cellidx_all,cellidx_all,pdcells);
        else
            cond_cellidx2 = cellidx_all(~ismember(cond_cellidx1,cellidx_all));
            condition{2} = ['non' condition{1}];
        end
        
        
        cellnr1 = length(cond_cellidx1); cellnr2 = length(cond_cellidx2);
        cellnr_pie(cellnr1,cellnr2,[],[1 1 0; [0.5 0.5 0.3]],condition);
        
    end
    
    
    
    set(H,'Position',get(0,'Screensize'));
    set(0, 'DefaultFigureRenderer', 'painters');
    
    
    % Save avg psth
    group_ev_dir2 = fullfile(group_ev_dir,'PSTH_avg'); fastif(~isdir(group_ev_dir2),mkdir(group_ev_dir2),0);
    
    if contains(condition{1},'bursting') || contains(condition{1},'rhythmic')
        
        fnm = fullfile(group_ev_dir2,[condition{1} num2str(cnr) '_' alignevent '_average.jpg']);   % jpg
        fnm2 = fullfile(group_ev_dir2,[condition{1} num2str(cnr) '_' alignevent '_average.fig']);   % fig
        fnm3 = fullfile(group_ev_dir2,[condition{1} num2str(cnr) '_' alignevent '_average.pdf']);   % fig
        
    else
        fnm = fullfile(group_ev_dir2,[alignevent '_average.jpg']);   % jpg
        fnm2 = fullfile(group_ev_dir2,[alignevent '_average.fig']);   % jpg
        fnm3 = fullfile(group_ev_dir2,[alignevent '_average.pdf']);  % fig
    end
    saveas(H,fnm);
    saveas(H,fnm2);
    
    close(H)
    
end













function cellnr_pie(cellnr1,cellnr2,cellnrRem,colors, labels)

% Plot pie chart with nr of cells

X = [cellnr1 cellnr2 cellnrRem];
P = pie(X);
patchHand = findobj(P, 'Type', 'Patch');

if ~isempty(cellnrRem)
    newColors = [colors; [0.6, 0.6, 0.6]];
    leglabels = {['n= ' num2str(X(1)) ' - ' labels{1}],...
        ['n= ' num2str(X(2)) '-' labels{2}],...
        ['n= ' num2str(X(3)) '-' labels{3}] };
else
    newColors = colors;
    leglabels = {['n= ' num2str(X(1)) ' - ' labels{1}],...
        ['n= ' num2str(X(2)) '-' labels{2}]};
    
end
set(patchHand, {'FaceColor'}, mat2cell(newColors, ones(size(newColors,1),1), 3))


lgd = legend(leglabels,'Location','bestoutside');
