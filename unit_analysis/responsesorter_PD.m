function [propname_stat, propname_boxstat] = responsesorter_PD(pdcells, resdir, tag, twin, issave)
%RESPONSESORTER_PD      Tests event-related activation/ inhibition 
%   [propname_stat, propname_boxstat] = responsesorter_PD(pdcells, resdir, tag, twin, issave)
%       -Plots PSTHs and raster plots in [-3 3] sec time window around a behavioral event of units listed in PDCELLS.
%       -Tests if there is a statistically significant activation/inhibition in test window (TWIN) compared 
%       to baseline period (see ULTIMATE_PSTH, BOXSTAT).
%       Results are saved also in CellBase (1 if activated; -1 if
%       inhibited; 0 if no significant change in firing rate).
% 
% Input parameters:
%     PDCELLS   cell array of unit IDs to analyse
%     
%     RESDIR	results directory to save PSTHs and raster plots
%     
%     TAG       character array, label of behavioral event
%     
%     TWIN      1x2 vector, test window in sec (ex. [0 1]) relative to event
%     
%     ISSAVE    1|0, if 1, save figures
% 
%  Output parameters:
%   PROPNAME_STAT       name of property in CellBase corresponding to ultimate
%                       psth stat testing (see ULTIMATE_PSTH)
%   
%   PROPNAME_BOXSTAT    name of property in CellBase corresponding to
%                       ranksum testing (see BOXSTAT)
%
% See also: ULTIMATE_PSTH, PARTITION_TRIALS, VIEWCELL2B, BOXSTAT

% Panna Hegedus, Hangya Balázs, Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


% Directories
if ~isfolder(resdir)
    mkdir(resdir);
end

numCells = length(pdcells);


% Load CellBase
load(getpref('cellbase','fname'),'CELLIDLIST');


% Add propoerty for grouping
pause(0.5)
propname_stat = [tag 'psth_stat_' num2str(twin)]
if ~ismember(propname_stat,listtag('prop'))
    insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname_stat)
end

propname_boxstat =[tag 'boxstat_' num2str(twin)];
if ~ismember(propname_boxstat,listtag('prop'))
    insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname_boxstat)
end

partition = 'all';   % partition trials
wn = [-4 4];   % full raster window in seconds
dt = 0.001;   % raster resolution, in seconds
sigma = 0.02;   % controls smoothing for 'spsth'
        
% Raster + PSTH
switch tag
    
    case 'StimulusOn'
        % Raster + PSTH aligned to stimulus onset
        alignevent = 'StimulusOn';   % trigger event
        shevent = 'KeyPress1';
        bwin = [-2.5 -1];   % baseline window for MW-test
        
        
    case 'KeyPress1'
        % Raster + PSTH aligned to stimulus onset
        alignevent = 'KeyPress1';   % trigger event
        shevent = 'StimulusOn';
        bwin = [-3 -1.5];   % baseline window for MW-test
        
        
    case 'StopSignal'
        % Raster + PSTH aligned to stimulus onset
        alignevent = 'StopSignal';   % trigger event
        shevent = 'StimulusOn';
        bwin = [-3 -1.5];   % baseline window for MW-test
        
    case 'Feedback'
        % Raster + PSTH aligned to stimulus onset
        alignevent = 'Feedback';   % trigger event
        shevent = 'StimulusOn';
        bwin = [1.5 3];   % baseline window for MW-test
end

% PSTH
for iC = 1:numCells
    cellid = pdcells{iC};   % current cell
    
    cellidt = regexprep(cellid,'\.','_');
    fnm = fullfile(resdir,[cellidt '_' alignevent '_' partition(2:end) '_boxstat.jpg']);
%     if exist(fnm)==2
%         fprintf('%s done.\n',cellid)
%         continue
%     else
%         fprintf('%s in progress.\n',cellid)
%     end
    
    try
        %%
        [binrast, stats] = rasterPSTH(cellid,alignevent,shevent,partition,wn,dt,sigma,bwin,twin,issave,resdir);
    catch
        fprintf('Error in rasterPSTH %s\n',cellid)
        continue
    end
    
    % Add property to CellBase - if both inhibited and activated, the
    % effect with larger change relative to baseline is considered

    if stats{1}.Wpi < 0.05 && stats{1}.Wpa < 0.05
        if diff([stats{1}.minvalue stats{1}.baseline])> diff([stats{1}.baseline stats{1}.maxvalue])
            st = setvalue(cellid,propname_stat,-1); 
        else
            st = setvalue(cellid,propname_stat,1); 
        end
    elseif stats{1}.Wpi < 0.05
            st = setvalue(cellid,propname_stat,-1);   % inhibited
    elseif stats{1}.Wpa < 0.05
            st = setvalue(cellid,propname_stat,1);   % activated
        else
            st = setvalue(cellid,propname_stat,0);   % non-responsive
        end

    if ~st
        error('MATLAB:vpresponsesorter:setvalueStatus','Error using setvalue.m')
    end
    
    % Find window indices
    time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector
    btime = bwin(1)*1000:dt*1000:bwin(2)*1000;   % baseline time vector
    baseline_inx = dsearchn(time',btime'); % baseline indices
    lbwin = abs(bwin(2)-bwin(1)); % length of baseline window in seconds
    
    ttime = twin(1)*1000:dt*1000:twin(2)*1000; % test time vector
    test_inx = dsearchn(time',ttime'); % test window indices
    ltwin = abs(twin(2)-twin(1)); % length of test window in seconds
    
    NumTrials = size(binrast,1);
    baseline_firing = nan(1,NumTrials);
    test_firing = nan(1, NumTrials);
    for t=1:NumTrials
        baseline_firing(t) = sum(binrast(t,baseline_inx))/lbwin;
        test_firing(t) = sum(binrast(t,test_inx))/ltwin;
    end
    [H, p]= boxstat(baseline_firing,test_firing, 'baseline', 'test', 0.05, 'paired')% statistics
    
    if (p < 0.05) && (median(baseline_firing)< median(test_firing)) % activation
        st = setvalue(cellid,propname_boxstat,1);
    elseif (p < 0.05) && (median(baseline_firing)> median(test_firing)) % inhibited
        st = setvalue(cellid,propname_boxstat,-1);
    elseif p > 0.05
        st = setvalue(cellid,propname_boxstat,0);   % non-responsive
    end
    
    
    if issave
        cellidt = regexprep(cellid,'\.','_');
        fnm = fullfile(resdir,[cellidt '_' alignevent '_' partition(2:end) '_boxstat.jpg']);
        saveas(H,fnm)
        fnm2 = fullfile(resdir,[cellidt '_' alignevent '_' partition(2:end) '_boxstat.fig']);
        saveas(H,fnm2)
        close all
    end
end

function [binrast, stats1] = rasterPSTH(cellid,alignevent,shevent,partition,wn,dt,sigma,bwin,twin,issave,response_resdir)

% Raster plot and PSTH
TE = loadcb(cellid,'TrialEvents');   % load trial events
SP = loadcb(cellid,'EVENTSPIKES');
fld = fieldnames(TE);
% if ~isequal(length(SP.event_stimes{1}),length(TE.(fld{1})))
%     error('MATLAB:vpresponsesorter:rasterPSTH:trialMismatch',...
%         'Trial number mismatch between TrialEvents and EVENTSPIKES.')
% end

% Raster plot
viewcell2b(cellid,'TriggerName',alignevent,'SortEvent','TrialStart','sigma',sigma,...
    'eventtype','behav','ShowEvents',{{shevent}},'Partitions',partition,'window',wn,'PSTHPlot',false,'dt',dt);
V_handle = gcf;
maximize_figure

% Peri-event time histogram
PSTHaxis_handle = findobj(allchild(V_handle),'type','axes','XLim',[0 1],'YLim',[0 1],'Tag','');   % handle for the empty PSTH axes
[~, ~, ~, ~, binrast, stats1] = ...
    ultimate_psth(cellid,'trial',alignevent,wn,...
    'dt',dt,'sigma',sigma,'parts',partition,'isadaptive',0,...
    'maxtrialno',Inf,'baselinewin',bwin,'testwin',twin,'relative_threshold',0.1,'display',true); % calculate psth

% % Plot & save
if ~iscell(stats1)
    stats1 = {stats1};   % only one PSTH
end
[~, tags] = partition_trials(TE,partition);
[mylabels, mycolors] = makeColorsLabels(@defineLabelsColors_Balazs,tags);
NumStats = length(stats1);
close_handle = nan(1,NumStats);
for iP =  1:NumStats
    Ls = findobj(allchild(stats1{iP}.axis_handle),'Type','line');   % all lines in the plot
    clr = findobj(Ls,'Type','line','Color','black');   % re-color the black PSTH
    set(clr,'Color',mycolors{iP})
    tx = findobj(allchild(stats1{iP}.axis_handle),'Type','text');   % re-position text object
    x_lim = xlim;
    y_lim = ylim;
    set(tx(1),'Position',[x_lim(1)+diff(x_lim)*0.1 (iP-1)*diff(y_lim)*0.4+y_lim(1)+diff(y_lim)*0.6 0])
    tx(1).String = regexprep(tx(1).String,'MW test',mylabels{iP});   % more meaningful labels
    if length(tx) > 1
        set(tx(2),'Position',[x_lim(1)+diff(x_lim)*0.1 (iP-1)*diff(y_lim)*0.4+y_lim(1)+diff(y_lim)*0.4  0])
        tx(2).String = regexprep(tx(2).String,'MW test',mylabels{iP});   % more meaningful labels
    end
    copyobj(tx,PSTHaxis_handle)   % copy PSTH and stat text to the raster figure
    copyobj(Ls,PSTHaxis_handle)
    xlabel(PSTHaxis_handle,['Time from ' alignevent])
    close_handle(iP) = stats1{iP}.figure_handle;
end

% Save figure
if issave
    cellidt = regexprep(cellid,'\.','_');
    fnm = fullfile(response_resdir,[cellidt '_' alignevent '_' partition(2:end) '.jpg']);
    saveas(V_handle,fnm)
    fnm2 = fullfile(response_resdir,[cellidt '_' alignevent '_' partition(2:end) '.fig']);
    saveas(V_handle,fnm2)
    fnm = fullfile(response_resdir,[cellidt '_' alignevent '_' partition(2:end) '.mat']);
    warning('off','MATLAB:Figure:FigureSavedToMATFile')
    save(fnm,'stats1')
end

close([V_handle, close_handle])