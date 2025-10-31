function respsort_partitions_PD(pdcells,alignevent,partition, twinds, resdir, issave,excl_beforeEv)
%RESPSORT_PARTITIONS_PD  Tests activation/ inhibition predicting the behavioral outcome
% RESPSORT_PARTITIONS_PD(pdcells,alignevent,partition, twinds, resdir, issave)
%       -Plots PSTHs and raster plots of partitioned trials in [-3 3] sec time window around 
%       a behavioral event of units listed in PDCELLS.
%       -Tests if there is a statistically significant activation/inhibition predicting 
%       behavioral outcome (see BOXSTAT) in test windows (TWINDS) by comparing median firing rate
%       of opposing partitions in predefined test windows (multiple test windows are defined). 
%       Results are saved also in CellBase:
%            1 if Partition nr 1. (Successful Stop)< Partition nr 2 (Failed Stop); 
%           -1 if Partition nr 1. (Successful Stop)> Partition nr 2 (Failed Stop); 
%            0 if no significant difference in firing rate
%
% Input parameters:
%     PDCELLS       cell array of unit IDs to analyse
%
%     ALIGNEVENT    character array, label of behavioral event
%
%     PARTITION     character array, label of partition criteria
%     
%     TWINDS        nx2 vector test windows in sec (ex. [-0.5 0; -1 0; -1.5 0; -2 0])
%     
%     RESDIR        results directory to save PSTHs and raster plots   
%     
%     ISSAVE        1|0, if 1, save figures
%
% See also: PULTIMATE_PSTH, PARTITION_TRIALS, VIEWCELL2B, BOXSTAT

% Panna Hegedus, Hangya Balázs, Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

dbstop if error

% Directories
if ~isfolder(resdir)
    mkdir(resdir);
end

numCells = length(pdcells);


% Load CellBase
load(getpref('cellbase','fname'),'CELLIDLIST');

if isempty(excl_beforeEv)
    resdirbox = fullfile(resdir,[alignevent '_' partition(2:end)], 'boxplots');
    resdirps = fullfile(resdir,[alignevent '_' partition(2:end)], 'raster_psths')
else
    resdirbox = fullfile(resdir,[alignevent '_' partition(2:end) '_ExclBef_' excl_beforeEv], 'boxplots');
    resdirps = fullfile(resdir,[alignevent '_' partition(2:end) '_ExclBef_' excl_beforeEv], 'raster_psths')
end
fastif(~isdir(resdirbox),mkdir(resdirbox),0);
fastif(~isdir(resdirps),mkdir(resdirps),0);


% Raster + PSTH aligned to stimulus onset

switch alignevent
    case 'StimulusOn'
        shevent = 'KeyPress1';
        wn = [-3 3];   % full raster window in seconds
        dt = 0.001;   % raster resolution, in seconds
        sigma = 0.02;   % controls smoothing for 'spsth'
        bwin = [-2 -1];   % baseline window for MW-test
    case 'StopSignal'
        shevent = 'StimulusOn';
        wn = [-3 3];   % full raster window in seconds
        dt = 0.001;   % raster resolution, in seconds
        sigma = 0.02;   % controls smoothing for 'spsth'
        if isempty(excl_beforeEv)
            bwin = [-3 -2];   % baseline window for MW-test
        else
            bwin = [-1.1 -0.1];
        end
        
    otherwise
        error('Add parameters for other events')
end


% PSTH

propname_stat = {};
for iC = 1:numCells
    
    
    
    cellid = pdcells{iC};   % current cell
    cellidt = regexprep(cellid,'\.','_');
    
    
    fprintf('%s...\n',cellid)
    
    % Load behav + spikes info
    TE = loadcb(cellid,'TrialEvents');   % load trial events
    SP = loadcb(cellid,'EVENTSPIKES');
    fld = fieldnames(TE);
    if ~isequal(length(SP.event_stimes{1}),length(TE.(fld{1})))
        error('MATLAB:vpresponsesorter:rasterPSTH:trialMismatch',...
            'Trial number mismatch between TrialEvents and EVENTSPIKES.')
    end
    
    
    if iC==1
        [~, tags] = partition_trials(TE,partition);
        [mylabels, mycolors] = makeColorsLabels(@defineLabelsColors_Balazs,tags);
    end
    
    % Peri-event time histogram
    [~, ~, ~, ~, binrast] = ...
        ultimate_psth(cellid,'trial',alignevent,wn,...
        'dt',dt,'sigma',sigma,'parts',partition,'isadaptive',0,...
        'maxtrialno',Inf,'baselinewin',bwin,'relative_threshold',0.1,'display',false,...
        'first_event',excl_beforeEv); % calculate psth
    
    
    viewcell2b(cellid,'TriggerName',alignevent,'SortEvent','TrialStart','sigma',sigma,...
        'eventtype','behav','ShowEvents',{{shevent}},'Partitions',partition,'window',wn,'PSTHPlot',true);
    psth_fig = gcf;
    set(psth_fig,'Visible','off')
    psth_ax = gca;
    
    % Compute stat form multiple windows
    sval = [];
    for tw = 1:size(twinds,1)
        
        
        twin = twinds(tw,:);
        
        
        % Add propoerty for grouping
        if iC==1
            propname_stat{tw} = ['x' alignevent '_' partition(2:end) '_stat_' num2str(twin)];
            
            if ~ismember(propname_stat{tw},listtag('prop'))
                insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname_stat{tw})
            end
            
        end
        
        
        
        % Find window indices
        time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector
        
        ttime = twin(1)*1000:dt*1000:twin(2)*1000; % test time vector
        test_inx = dsearchn(time',ttime'); % test window indices
        ltwin = diff(twin); % length of test window in seconds

        
        if ~iscell(binrast)
            binrast = {binrast};
        end
        
        
        if length(binrast)<2
            continue
        end
        
        NumTrials1 = size(binrast{1},1);
        NumTrials2 = size(binrast{2},1);
        
        part1_firing = nan(1,NumTrials1);
        part2_firing = nan(1, NumTrials2);
        
        for t=1:NumTrials1
            part1_firing(t) = sum(binrast{1}(t,test_inx))/ltwin; % StopPartition =1 = FailedStop
        end
        
        for t=1:NumTrials2
            part2_firing(t) = sum(binrast{2}(t,test_inx))/ltwin; % StopPartition =2 = SuccessfulStop
        end
        
        
        [H, p]= boxstat(part1_firing,part2_firing,mylabels{1}, mylabels{2}, 0.05, 'nonpaired')% statistics
       
        
        
        pos = get(psth_ax,'Position');
        xL = get(psth_ax,'XLim');
        annotXpos= interp1([xL(1) xL(end)],[pos(1) pos(1)+pos(3)],twin);
        annotYpos = [pos(2)+0.1+tw*0.025 pos(2)+0.1+tw*0.025];
        
        if p<0.05
            annotcol = 'red';
        else
            annotcol = 'black';
        end
        
        annotation(psth_fig,'line',annotXpos,annotYpos,'Color',annotcol)
        annotation(psth_fig,'textbox',[annotXpos(1) annotYpos(1)-0.025 annotXpos(2) annotYpos(1)],'String',['p = ' num2str(p)],...
            'Color',annotcol,'LineStyle','none','VerticalAlignment','bottom')
        
        
        
        if (p < 0.05) && (median(part1_firing)< median(part2_firing)) % Failed < Succ   
            
            st = setvalue(cellid,propname_stat{tw},-1);
        elseif (p < 0.05) && (median(part1_firing)> median(part2_firing)) % Failed > Succ    
            
            st = setvalue(cellid,propname_stat{tw},1);
        else
            
            st = setvalue(cellid,propname_stat{tw},0);  % non-responsive
        end
        
        if issave
            fnm = fullfile(resdirbox,[cellidt '_' alignevent '_' partition(2:end) '_boxstat_' num2str(twin)]);
            if ~isempty(excl_beforeEv)
                fnm = [fnm '_ExclBef_' excl_beforeEv];
            end
            saveas(H,[fnm '.jpg'])
            saveas(H,[fnm '.fig'])
            close(H)
        end
        
    end
    
    set(psth_fig,'Position',get(0,'Screensize'));
    
    
    if issave
        
        fnmp = fullfile(resdirps,[cellidt '_' alignevent '_' partition(2:end) '_psth']);
        if ~isempty(excl_beforeEv)
            fnmp = [fnmp '_ExclBef_' excl_beforeEv];
        end
        saveas(psth_fig,[fnmp '.fig'])
        saveas(psth_fig,[fnmp '.jpg'])
        
        close(psth_fig)
    end
end



