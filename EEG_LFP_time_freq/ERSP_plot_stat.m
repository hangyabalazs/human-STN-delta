function ERSP_plot_stat(sess2analyse,EventTypes,SubEventTypes, epoch_win,baseline_win,varargin)
%ERSP_PLOT_STAT     Event-related spectral perturbation (ERSP) maps with statistics
%   ERSP_plot_stat(sess2analyse,EventTypes,SubEventTypes,...)
%       - plots ERSP maps patient-by-patient using normalized epoched data (see time_freq_patients.m)
%       - calculates time-frequency patient averages
%       - plots ERSP maps of patient averages
%       - if postoperative EEG, draws topoplots using averages of across specified time-frequency windows
%
%
% Required inputs:
%     SESS2ANALYSE      struct containing all necessary information (name of patient, side
%                       of experiment, tag of condition, session folder path) of
%                       session data that need to be analysed (see getdata2analyse)
%
%     EVENTTYPES        1xN cell array of event labels, ex: {'StimulusOn','StopSignal','KeyPress1','Feedback'};
%
%     SUBEVENTTYPES     Nx2 cell array of partition ("subevent") labels, each row
%                       corresponds to an event label, columns to partitions
%                       {'CueStim','StopStim';'FailedStopTrial','SuccesfulStopTrial';'CueResponse','StopResponse';'Correct','Error';};
%
%     EPOCH_WIN         1x2 vector, time window relative to event timestamp in sec, for data epoching (ex: [-2 2])
%
%     BASELINE_WIN       1x2 vector, time window relative to event timestamp in sec, for baseline correction
%
% Optional inputs:
%   'stat_time'     1x2 vector, time window relative to event timestamp in
%                   sec, to perform statistics (default value: [-1 1])
%
%   'csd'           true | false, if true, CSD transformed EEG data is used
%                   (relevant only for postoperative EEG data) (default value: true)
%
%   'bipol'         true | false, if true, F4-F3 bipolar derivation of EEG data is used
%                   (relevant only for postoperative EEG data) (default value: false)
%
%   'alphas'        numeric or vector, list of alpha threshold to consider when performing statistics
%                   ERSP masks (1 where stat. significant, 0 otherwise) are created using the the specified alpha threshold
%                   if NaN, no statistical test is performed
%                   (default value: NaN)
%
%   'patgroup_nm'   cell array of patient group labels
%               {}, includes all patients
%               {'tremor-dominant','akinetic-rigid','mixed'}, draws separate plots for each clinical group
%               {'RTdecrease','RTincrease'}, - draws separate plots for each patient group based on preop-postop RT change
%               (default value: {})
%
%   'bypatient'     if true, ERSP maps are generated for each patient
%                   true | false (default value: true)
%
%   'patientavg'    if true, ERSP maps are generated for patient averages
%                   true | false (default value: true)
%
%   'ersp_freq'     1x2 vector, frequency limits of ERSP map (default value: [1 80])
%
%   'subregion'     character array or cell array, plots channel data
%                   derived from listed STN subregion
%                   (relevant only for intraop LFP data)
%          'all' | 'Motor' | 'Limbic' | 'Associative'  (default value: 'all')
%
%   'close2centr'   character array or cell array, plots channel data
%                   closest to listed STN subregion
%                   (relevant only for intraop LFP data)
%          'all' | 'Motor' | 'Limbic' | 'Associative'  (default value: 'all')
%
%   'chanmean'      if true, channel averaged LFP data is used, if false
%                   LFP data is plotted channel-by-channel
%                   (relevant only for intraop LFP data)
%                   1 | 0 (default value: 1)
%
%   'topo_freq_nms' cell array, list of frequency band labels to represent
%                   on topoplots (relevant only for postop EEG)
%                   (default value: {'delta'})
%
%   'topot_freqs'   matrix, frequency band limits, each row corresponds to a label
%                   in 'topo_freq_nms' (relevant only for postop EEG)
%                   (power averaged across frequency components within limits)
%
%   'topobin'       integer, time window represented on one topoplot in ms
%                   (power averaged within time bin) (relevant only for postop EEG)
%
% See also: TIME_FREQ_PATIENTS

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

prs = inputParser;
addRequired(prs,'sess2analyse',@isstruct);
addRequired(prs,'EventTypes',@iscell);
addRequired(prs,'SubEventTypes',@iscell);
addRequired(prs,'epoch_win',@isvector);
addRequired(prs,'baseline_win',@isvector);
addParameter(prs,'side','left',@(x) ischar(x)||iscell(x));
addParameter(prs,'condition','stimoff',@(x) ischar(x)||iscell(x));
addParameter(prs,'stat_time',[-1 1],@isvector);
addParameter(prs,'csd',true,@islogical);
addParameter(prs,'bipol',false,@islogical);
addParameter(prs,'alphas',NaN,@(x) isnumeric(x)||isvector(x));
addParameter(prs,'patgroup_nm',{},@iscell);
addParameter(prs,'bypatient',true,@islogical);
addParameter(prs,'patientavg',true,@islogical);
addParameter(prs,'ersp_freq',[1 80],@isvector);
addParameter(prs,'subregion','all',@(x) ischar(x)||iscell(x));
addParameter(prs,'close2centr','all',@(x) ischar(x)||iscell(x));
addParameter(prs,'chanmean',1,@isnumeric);
addParameter(prs,'topo_freq_nms',{'delta'},@iscell);
addParameter(prs,'topo_freqs',[1 4],@ismatrix);
addParameter(prs,'topobin',500,@isnumeric);
parse(prs,sess2analyse,EventTypes,SubEventTypes,epoch_win,baseline_win,varargin{:})
p = prs.Results;


rectype = sess2analyse(1).rectype;

if strcmp(rectype, 'LFP')
    p.csd = false;
    p.bipol = false;
end


    
    
    % Power averages across patients
    save_avgepoch_pow2(sess2analyse,'EventTypes',EventTypes,'SubEventTypes',SubEventTypes,...
        'subregion',p.subregion,'close2centr',p.close2centr,'csd',p.csd,'bipol',p.bipol,'measure','Pow','alpha',p.alphas,...
        'patgroup_nm',p.patgroup_nm, 'baseline_win',baseline_win,'stat_win',p.stat_time,...
        'side',p.side,'condition',p.condition,'chanmean',p.chanmean);
    
    
    gr_label = {};
    if ~isempty(p.patgroup_nm)
        for k = 1:length(p.patgroup_nm)
            gr_label{k} = [p.side '_' p.condition '_RACEfalse_'  p.patgroup_nm{k}];
        end
    else
        gr_label{1} = [p.side '_' p.condition '_RACEfalse'];
    end
    
    
    % Plot ERSP patient-by-patient & patient AVG
    PD_ERSPs(sess2analyse,epoch_win,'EventTypes',EventTypes,'SubEventTypes',SubEventTypes,'Subregion',p.subregion,...
        'bypatient',p.bypatient,'patientavg',p.patientavg,'freq_range',p.ersp_freq,'time_range',p.stat_time,'chanmean',p.chanmean,...
        'grcond_names',gr_label,'stat',true,'measure','Pow','csd',p.csd,'bipol',p.bipol,'baseline_win',baseline_win,'cLim',[-1 1]);
    
    
    
    % Topoplots for each patient (event-averages)
    ssrt_topoplot(sess2analyse,'plot_win',p.stat_time,'topo_freq_nms',p.topo_freq_nms,'topo_freqs',p.topo_freqs,...
        'topobin',p.topobin,'cLim',[-2 2],'individual',p.bypatient,'patavg',p.patientavg,'csd',p.csd,'measure','Pow',...
        'stat',true,'grcond_names',gr_label,'baseline_win',baseline_win,'EventTypes',EventTypes,'SubEventTypes',SubEventTypes);
    
    
   
end



%--------------------------------------------------------------------------
function PD_ERSPs(sess2analyse,epoch_win, varargin)

% Saves average time-freq maps into AvgPow matrix (results directory: figdir\TF)
global rootdir figdir_pd

prs = inputParser;
addRequired(prs,'sess2analyse',@isstruct);
addRequired(prs,'epoch_win',@isvector);
addParameter(prs,'rectype','',@ischar);
addParameter(prs,'EventTypes',{'StimulusOn','StopSignal'},@iscell);
addParameter(prs,'SubEventTypes', {'FailedStopTrial','SuccesfulStopTrial';...
    'FailedStopTrial','SuccesfulStopTrial'},@iscell);
addParameter(prs,'preproc','epochs',@ischar);
addParameter(prs,'chanmean',1,@isnumeric);
addParameter(prs,'subregion','all',@ischar);
addParameter(prs,'close2centr','all',@ischar);
addParameter(prs,'maxpower','all',@ischar);
addParameter(prs,'bypatient',true,@islogical);
addParameter(prs,'csd',false,@islogical);
addParameter(prs,'bipol',false,@islogical);
addParameter(prs,'patientavg',true,@islogical);
addParameter(prs,'freq_range',[0.5 100],@isvector);
addParameter(prs,'time_range',[-2 2],@isvector);
addParameter(prs,'baseline_win',[],@(x) isvector(x)|isnumeric(x));
addParameter(prs,'grcond_names',{'left_stimoff','right_stimoff'},@iscell); % {'left_stimoff','left_stimon','right_stimoff','right_stimon'}
addParameter(prs,'stat',false,@islogical);
addParameter(prs,'measure','Pow',@ischar);
addParameter(prs,'cLim',[],@(x) isnumeric(x)|isvector(x));


parse(prs,sess2analyse,epoch_win,varargin{:});
p = prs.Results;


%%
% Recording type
if isempty(p.rectype)
    rectype = sess2analyse(1).rectype;
else
    rectype = p.rectype;
end

% Recordint time
rectime = sess2analyse(1).rectime;
if ~strcmp(rectime,'postop')
    p.csd = false; p.bipol = false;
end

% Load subregion data if necessary
if ~strcmp(p.subregion,'all') && strcmp(rectype,'LFP')
    load(fullfile(rootdir,'STN_loc','STN_loc.mat'));
    subreg_names = STN_loc.subreg_names;
    sinx = find(strcmp(subreg_names,p.subregion));
end

% Load freq components
load(fullfile(rootdir,'freq_components.mat'))
f_ind = intersect(find(f>=p.freq_range(1)),find(f<=p.freq_range(2)));

% Time
new_sr = 50;

if ~isequal(p.time_range,p.epoch_win)
    times = p.epoch_win(1):1/new_sr:p.epoch_win(2);
    new_times = p.time_range(1):1/new_sr:p.time_range(2);
    newtinx = dsearchn(times',new_times'); newtinx = newtinx(1:end-1);
else
    newtinx = 1:200;
end

if ~isempty(p.baseline_win)
    
    baslims = dsearchn(times',p.baseline_win');
    basnm = ['BASLims' num2str(baslims(1)) '_' num2str(baslims(2))]; basnm2 = basnm;
    
    if newtinx(1)>=baslims(1)
        newtinx = newtinx-baslims(1)+1; % bas-corrected data is saved from start of baseline until end of original data epoch
    else
        fprintf('Start of time period before baseline start\n');
        return;
    end
    timerang = [max( [ p.baseline_win(1) p.time_range(1) ]) p.time_range(2)];
else
    basnm = ['TRN_AVG'];
    basnm2 = 'TrialNorm';
    timerang = p.time_range;
end


if isempty(p.SubEventTypes); subnr = 3; else; subnr = 1:3; end;



%% Average over epochs -patient by patient
if p.bypatient
    for snr = 1:length(sess2analyse)
        curr_resdir = sess2analyse(snr).folder;
        side = sess2analyse(snr).side;
        tag = sess2analyse(snr).tag;
        patnm = sess2analyse(snr).patient;
        
        
        fprintf('%s %s %s...\n', patnm,side, tag);
        
        % Load event averages for patients
        if p.bipol
            resdir = fullfile(curr_resdir, 'EventAVGs_bipol');
        else
            resdir = fullfile(curr_resdir, ['EventAVGs_CSD' char(string(p.csd))]);
        end
        
        
        
        
        % Channels
        if strcmp(rectype,'LFP') && p.chanmean==1
            channels = {'chanmean'};
        elseif strcmp(rectype,'EEG') && strcmp(rectime,'intraop')
            channels = {'F4'}; p.chanmean = 0;
        elseif strcmp(rectype,'EEG') && strcmp(rectime,'postop')
            if ~p.bipol
                load(fullfile(rootdir,'postop_EEG_chanlocs.mat'))
                channels = {chanlocs.labels};
            else
                channels = {'F4-F3'};
            end
        end
        
        
        
        % Loop over events
        for ei = 1:length(p.EventTypes)
            event = p.EventTypes{ei};
            
            
            
            tf = [];
            % Loop over subevents
            for subei = subnr
                if subei<3
                    subevent = p.SubEventTypes{ei,subei};
                else
                    subevent = event;
                end
                
                
                if strcmp(rectype,'LFP') && p.chanmean==1
                    
                    evavg_nm = [ event '_' subevent '_chan_AVGs.mat'];
                    evSTAT_nm = [ event '_' subevent '_chan_STATs.mat'];
                else
                    evavg_nm = [ event '_' subevent '_AVGs.mat'];
                    evSTAT_nm = [ event '_' subevent '_STATs.mat'];
                end
                
                
                try
                    load(fullfile(resdir,evavg_nm))
                    load(fullfile(resdir,evSTAT_nm))
                catch
                    fprintf('No files %s %s %s %s %s \n', patnm, side, tag, event, subevent)
                    continue;
                end
                
                
                
                if strcmp(rectype,'LFP') && p.chanmean~=1
                    channels = fieldnames(EventAVGs);
                end
                
                % Loop over channels
                for ch = 1:length(channels)
                    act_chan = channels{ch};
                    if p.bipol
                        act_chan(strfind(act_chan, '-')) = '_';
                    end
                    
                    % Result directory
                    if ~p.bipol
                        resdir_fig = fullfile(figdir_pd,[rectime '_' rectype], [p.measure '_ERSP'],['ERSP_patients_CSD' char(string(p.csd))],[event '_' subevent]);
                    else
                        resdir_fig = fullfile(figdir_pd,[rectime '_' rectype], [p.measure '_ERSP'],'ERSP_patients_bipol',[event '_' subevent]);
                    end
                    
                    if ~isdir(resdir_fig); mkdir(resdir_fig); end;
                    
                    
                    tf = EventAVGs.(act_chan).(basnm);
                    
                    
                    if isempty(tf)
                        fprintf('No %s AT ALL: %s %s\n',p.subregion,patnm,side)
                        continue
                    end
                    
                    
                    % FIGURE
                    P = tf(f_ind,newtinx);
                    fig = figure;
                    [~,~,~] = spectr_fig(P,f(f_ind),new_sr,p.freq_range(2),p.freq_range(1),timerang, [],0, p.cLim);
                    
                    load(fullfile(resdir,evSTAT_nm))
                    if p.stat
                        mask_ersp = EventSTAT.(act_chan).(basnm2).mask_ersp;
                        hold on;
                        contour(mask_ersp(f_ind,newtinx),'Color','white');
                    end
                    
                    
                    % Title
                    try
                        trnr = EventSTAT.(act_chan).(basnm2).trialnr;
                    catch
                        trnr = 0;
                    end
                    title({['Event:',subevent,',',' Channel:' act_chan],[patnm ', ' side ', ' tag],['n=' num2str(trnr) ' trials']});
                    
%                     fig = gcf;
                    
                    % Save figure
                    fnm = [resdir_fig filesep patnm '_' side '_' tag ...
                        '_' act_chan '_' p.subregion '_FREQ' num2str(p.freq_range) '_WIN' num2str(p.time_range) '_BAS' num2str(p.baseline_win)];
                    saveas(fig, [ fnm '.png']);
                    saveas(fig,[ fnm  '.fig']);
                    close(fig)
                    
                    
                    
                    
                    
                    
                end
            end
        end
    end
end
%% Average over patients

if p.patientavg
    % Channel names
    if strcmp(rectype,'EEG') && strcmp(rectime,'postop')
        if ~p.bipol
            load(fullfile(rootdir,'postop_EEG_chanlocs.mat')); choi = {chanlocs.labels}; clear chanlocs;
        else
            choi = {'F4-F3'};
        end
    elseif strcmp(rectype,'EEG') && strcmp(rectime,'intraop')
        choi = {'F4'};
        
        
    elseif strcmp(rectype,'LFP') % choose channels for averagin...
        if ~strcmp('all',p.subregion) % ...channels from spec. subregion
            choi = {['chanmean-' p.subregion]};
        elseif ~strcmp('all',p.close2centr) % ...closest channels to spec. subregion centroid
            choi = {['chans_close2_' p.close2centr]};
        elseif ~strcmp('all',p.maxpower) % ...channels with max power within spec. freq. band
            choi = {['chans-max-' p.maxpower]};
        else
            choi = {'chanmean'}; % ... all channels
        end
    end
    
    
    
    if p.bipol
        resdir =fullfile(figdir_pd,[rectime '_' rectype],'TF_bipol');
    elseif strcmp(rectime, 'postop')
        resdir =fullfile(figdir_pd,[rectime '_' rectype],['TF_CSD' char(string(p.csd))]);
    else
        resdir =fullfile(figdir_pd,[rectime '_' rectype],'TF');
    end
    if strcmp(p.measure,'Pow');
        avgnm = ['avgepoch_CSD' char(string(p.csd))];
    end;
    
    
    
    % Loop over groups (left_stimoff, left_stimon etc.)
    for gri = 1:length(p.grcond_names)
        grnm = p.grcond_names{gri};
        
        % Loop over events
        for ei = 1:length(p.EventTypes)
            event = p.EventTypes{ei};
            
            % Loop over subevents
            for sei = subnr
                if sei<3
                    subevent = p.SubEventTypes{ei,sei};
                else
                    subevent = event;
                end
                
                % Result directory
                if ~p.bipol
                    resdir_fig = fullfile(figdir_pd,[rectime '_' rectype], [p.measure '_ERSP'],['ERSP_patientAVG_CSD' char(string(p.csd))],grnm,[event '_' subevent]);
                else
                    resdir_fig = fullfile(figdir_pd,[rectime '_' rectype], [p.measure '_ERSP'],'ERSP_patientAVG_bipol',grnm,[event '_' subevent]);
                    
                end
                
                if ~isdir(resdir_fig); mkdir(resdir_fig); end;
                
                close all
                
                % Load average matrix by channel/ group of channels
                for ch = 1:length(choi)
                    act_chan = choi{ch};
                    if p.bipol
                        act_chan(strfind(act_chan, '-')) = '_';
                    end
                    resdir_chan = fullfile(resdir,act_chan); if ~isdir(resdir_chan); mkdir(resdir_chan); end;
                    
                    
                    
                    avgepoch_pow = [];
                    load(fullfile(resdir_chan,[event '_' subevent '_' avgnm '_pow_' grnm  num2str(p.baseline_win) '.mat']));
                    
                    
                    % FIGURE
                    [~,~,~] = spectr_fig(avgepoch_pow(f_ind,newtinx),f(f_ind),new_sr,p.freq_range(2),p.freq_range(1),p.time_range,[], 0,p.cLim);
                    
                    
                    load(fullfile(resdir_chan,[event '_' subevent '_' avgnm '_STAT_' grnm  num2str(p.baseline_win) '.mat']));
                    if p.stat
                        mask_ersp = avgepoch_STAT.mask_ersp;
                        hold on;
                        contour(mask_ersp(f_ind,newtinx),'Color','white');
                    end
                    
                    % Title
                    patnr = avgepoch_STAT.patientnr;
                    
                    title({['Event:',subevent,',',' Channel:' act_chan],['n=' num2str(patnr) ' patients'],...
                        ['Baseline: ' num2str(p.baseline_win)]});
                    
                    fig = gcf;
                    
                    % Save figure
                    fnm = [resdir_fig filesep 'PatientAVG_' act_chan '_' p.subregion '_' num2str(p.freq_range)  '_' num2str(p.time_range) '_' num2str(p.baseline_win)];
                    saveas(fig, [ fnm '.png']);
                    saveas(fig,[ fnm '.fig']);
                    saveas(fig,[ fnm '.pdf']);
                    close(fig)
                    
                    
                end
            end
        end
    end
end
end





%--------------------------------------------------------------------------
function save_avgepoch_pow2(sess2analyse, varargin)

% Saves average time-freq maps into AvgPow matrix (results directory: figdir\TF)
global rootdir figdir_pd

prs = inputParser;
addRequired(prs,'sess2analyse',@isstruct);
addParameter(prs,'rectype','',@ischar);
addParameter(prs,'EventTypes',{'StimulusOn','StopSignal'},@iscell);
addParameter(prs,'SubEventTypes', {'FailedStopTrial','SuccesfulStopTrial';...
    'FailedStopTrial','SuccesfulStopTrial'},@iscell);
addParameter(prs,'side','left',@ischar);
addParameter(prs,'condition','stimoff',@ischar);
addParameter(prs,'patgroup_nm',{},@iscell);
addParameter(prs,'goodrace',false,@islogical);

addParameter(prs,'csd',false,@islogical);
addParameter(prs,'bipol',false,@islogical);
% addParameter(prs,'preproc','epochs',@ischar);
addParameter(prs,'chanmean',1,@isnumeric);
addParameter(prs,'subregion','all',@(x) ischar(x)||iscell(x)); % 'all' | 'no' | 'Motor' | ...
addParameter(prs,'close2centr','all',@ischar); % 'all' | 'Motor' | ...
addParameter(prs,'maxpower','all',@ischar); % 'all' | 'beta'
addParameter(prs,'iscell',[],@(x) isvector(x)||islogical(x));
addParameter(prs,'measure','Pow',@ischar);
addParameter(prs,'alpha',NaN,@isnumeric);
addParameter(prs,'mcorrect','fdr',@ischar); % 'fdr' | 'none'
addParameter(prs,'baseline_win',[],@(x) isvector(x)||isnumeric(x));
addParameter(prs,'stat_win',[],@(x) isvector(x)||islogical(x));




parse(prs,sess2analyse,varargin{:});
p = prs.Results;


%% Directories & filenames

rectype = sess2analyse(1).rectype;
rectime = sess2analyse(1).rectime;

if ~strcmp(rectime,'postop')
    p.csd = false; p.bipol = false;
end

if strcmp(p.measure,'Pow')
    
    if p.bipol
        tfnm = 'TF_bipol';
    elseif strcmp(rectime,'postop')
        tfnm = ['TF_CSD' char(string(p.csd))];
    else
        tfnm = 'TF';
    end
    
    res_nm = ['avgepoch_CSD' char(string(p.csd))];
    
    % elseif strcmp(p.measure,'Phas')
    %     if ~p.csd; tfnm = 'TF';  evavg_nm = 'ITPC.mat';   res_nm = 'avgITPC';
    %     else; tfnm = 'TF_CSD';  evavg_nm = 'ITPC_CSD.mat';   res_nm = 'avgITPC_CSD';
    %     end;
end


resdir =fullfile(figdir_pd,[rectime '_' rectype],tfnm);
if ~isdir(resdir); mkdir(resdir); end;


if strcmp(rectype,'EEG') && strcmp(rectime,'postop')
    if ~p.bipol
        load(fullfile(rootdir,'postop_EEG_chanlocs.mat')); choi = {chanlocs.labels}; clear chanlocs;
    else
        choi = {'F4-F3'};
    end
elseif strcmp(rectype,'EEG') && strcmp(rectime,'intraop')
    choi = {'F4'};
    
    
elseif strcmp(rectype,'LFP') % choose channels for averaging: ...
    if ~strcmp('all',p.subregion) % ...channels from spec. subregion
        choi = {['chanmean-' p.subregion]};
    elseif ~strcmp('all',p.close2centr) % ...closest channels to spec. subregion centroid
        choi = {['chans_close2_' p.close2centr]};
    elseif ~strcmp('all',p.maxpower) % ...channels with max power within spec. freq. band
        choi = {['chans-max-' p.maxpower]};
    else
        choi = {['chanmean']}; % ... all channels
    end
end

%% Select patient group

% Select sessions according to side
sidinx = ismember({sess2analyse.side},p.side);
s2a_new = sess2analyse(sidinx);


% Select sessions according to condition
coinx = ismember({s2a_new.tag},p.condition);
s2a_new = s2a_new(coinx);


% Sort sessions according to clinical group
s2a_groups = {}; group_lab = {};
if ~isempty(p.patgroup_nm)
    patgroups = clinical_groups(p.patgroup_nm);
    for k = 1:length(patgroups)
        patinx = ismember({s2a_new.patient},patgroups{k});
        s2a_groups{k} = s2a_new(patinx);
        group_lab{k} = [p.side '_' p.condition '_RACE' char(string(p.goodrace)) '_'  p.patgroup_nm{k}];
    end
else
    s2a_groups{1} = s2a_new;
    group_lab{1} = [p.side '_' p.condition '_RACE' char(string(p.goodrace))];
end


times = -2:1/50:2; times = times(1:end-1);

if ~isempty(p.baseline_win)
    
    baslims = dsearchn(times',p.baseline_win');
    basinx = baslims(1):baslims(2);
    %     statlims = dsearchn(times',p.stat_win');
    %     statlims = statlims-baslims(1)+1;
    %     basnm = ['BASLims' num2str(baslims(1)) '_' num2str(baslims(2))];
else
    basinx = 1:length(times);
end


load(fullfile(rootdir,'freq_components.mat'));
%%
if isempty(p.SubEventTypes); subnr = 3; else; subnr = 1:3; end;

for k = 1:length(s2a_groups)
    gr_label = group_lab{k};
    
    
    % Collect averages for each patient
    AP = collect_avgs(s2a_groups{k},p.csd,p.bipol,p.EventTypes,p.SubEventTypes,choi,p.subregion,p.close2centr,p.maxpower,p.chanmean,'TRN_AVG');
    
    for ci = 1:length(choi)
        
        act_chan = choi{ci};
        
        if contains(act_chan, '-')
            act_chan(strfind(act_chan,'-')) = '_' ;
        end
        
        AvgPow = AP.(act_chan).AvgPow;
        
        resdir_chan = fullfile(resdir,act_chan); if ~isdir(resdir_chan); mkdir(resdir_chan); end;
        for eii = 1:length(p.EventTypes)
            
            event = p.EventTypes{eii};
            
            for ii = subnr
                if ii<3
                    subevent = p.SubEventTypes{eii,ii};
                else
                    subevent = event;
                end
                
                
                
                
                evg_epochs = cat(3,AvgPow.(event).(subevent){:});
                
                
                %               evg_epochs2 = evg_epochs(:,statlims(1):statlims(2),:);
                
                if ~isempty(p.baseline_win)
                    evg_epochs2 = evg_epochs(:,baslims(1):end,:);
                    
                    B2 = repmat(nanmean(evg_epochs(:,basinx,:),[2 3]),[1 size(evg_epochs2,2)]);
                    Bsd2 = repmat(std(evg_epochs(:,basinx,:),[],[2 3],'omitnan'),[1 size(evg_epochs2,2)]);
                    
                    
                    
                    B_AVG2 = (nanmean(evg_epochs2,3)- B2)./ Bsd2; % Baseline correction for trial average for visualization
                    avgepoch_pow = B_AVG2;
                else
                    evg_epochs2 = evg_epochs;
                    avgepoch_pow = mean(evg_epochs,3);
                end
                
                
                % Stat
                
                if ~isnan(p.alpha)
                    hold on;
                    [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(evg_epochs2,f,p.alpha,1000,false,p.mcorrect,[],1:length(basinx));
                    
                    avgepoch_STAT.p_ersp =  exactp_ersp;
                    avgepoch_STAT.mask_ersp =  maskersp;
                    avgepoch_STAT.alphafdr =  alphafdr;
                    avgepoch_STAT.alpha =  p.alpha;
                    avgepoch_STAT.patientnr =  size(evg_epochs,3);
                end
                
                
                
                
                
                save(fullfile(resdir_chan,[event '_' subevent '_' res_nm '_pow_' gr_label num2str(p.baseline_win) '.mat']),'avgepoch_pow');
                save(fullfile(resdir_chan,[event '_' subevent '_' res_nm '_STAT_' gr_label num2str(p.baseline_win) '.mat']),'avgepoch_STAT');
                avgepoch_pow = []; avgepoch_STAT = struct;
            end
            
        end
    end
end
end





%--------------------------------------------------------------------------
function AP = collect_avgs(sess2analyse,csd,bipol,EventTypes,SubEventTypes,choi,subregion,close2centr,maxpower,chanmean,basnm)

global rootdir
if isempty(SubEventTypes); subnr = 3; else; subnr = 1:3; end;
rectype = sess2analyse(1).rectype;


if ~strcmp(subregion,'all')
    load(fullfile(rootdir,'STN_loc','STN_loc.mat'));
    subregs = STN_loc.subreg_names;
    sinx = ismember(lower(subregs),lower(subregion));
end

AP = struct;
for snr = 1:length(sess2analyse)
    curr_resdir = sess2analyse(snr).folder;
    side = sess2analyse(snr).side;
    tag = sess2analyse(snr).tag;
    patnm = sess2analyse(snr).patient;
    
    
    fprintf('%s %s %s...\n', patnm,side, tag);
    
    if bipol
        
        resdir = fullfile(curr_resdir, 'EventAVGs_bipol');
    else
        resdir = fullfile(curr_resdir, ['EventAVGs_CSD' char(string(csd))]);
    end
    
    for ei = 1:length(EventTypes)
        event = EventTypes{ei};
        
        
        for subei = subnr
            if subei<3
                subevent = SubEventTypes{ei,subei};
            else
                subevent = event;
            end
            
            if strcmp(rectype,'LFP') && chanmean==1
                evavg_nm = [ event '_' subevent '_chan_AVGs.mat'];
            else
                evavg_nm = [ event '_' subevent '_AVGs.mat'];
            end
            
            try
                load(fullfile(resdir,evavg_nm));
                
            catch
                fprintf('No files %s %s %s.\n', patnm,side, tag);
                continue
            end
            
            
            for ci = 1:length(choi)
                
                act_chan = choi{ci};
                
                if contains(act_chan, '-')
                    act_chan(strfind(act_chan,'-')) = '_' ;
                end
                
                if ~isfield(AP,act_chan)
                    AP.(act_chan).AvgPow = [];
                end
                
                if ~isfield(AP.(act_chan).AvgPow,event) % create new side-tag group field
                    AP.(act_chan).AvgPow.(event) = [];
                end
                
                
                if ~isfield(AP.(act_chan).AvgPow.(event),subevent); AP.(act_chan).AvgPow.(event).(subevent) = {}; end;
                
                if ~strcmp(rectype,'LFP') || (strcmp(rectype,'LFP') && chanmean==1)
                    
                    
                    
                    if isempty(AP.(act_chan).AvgPow.(event).(subevent))
                        AP.(act_chan).AvgPow.(event).(subevent){1} = EventAVGs.(act_chan).(basnm);
                    else
                        AP.(act_chan).AvgPow.(event).(subevent){end+1} = EventAVGs.(act_chan).(basnm);
                    end
                    
                    
                    
                    
                    
                elseif strcmp(rectype,'LFP')
                    
                    lfpchans = fieldnames(EventAVGs);
                    
                    lfpavg = {};
                    for lfpci = 1:length(lfpchans)
                        lfp_act_chan = lfpchans{lfpci};
                        
                        
                        
                        if isempty(EventAVGs.(lfp_act_chan).(basnm))
                            fprintf('No %s channel.\n',lfp_act_chan);
                            continue
                        end
                        
                        
                        if contains(lower(subregion),'no')
                            if any(~ismember(STN_loc.Patients.(patnm).(side).Channels.(lfp_act_chan).Subregion,-1));
                                fprintf('Not %s: %s %s %s\n',subregion,patnm,side,lfp_act_chan)
                                continue
                            end
                        elseif ~strcmp(subregion,'all')
                            if ~any(ismember(STN_loc.Patients.(patnm).(side).Channels.(lfp_act_chan).Subregion(sinx),[1 2]));
                                fprintf('Not %s: %s %s %s\n',subregion,patnm,side,lfp_act_chan)
                                continue
                            end
                        elseif ~strcmp('all',close2centr)
                            load(fullfile(rootdir,'Channels_close2centroids.mat'));
                            mychan = centr_tab{[patnm '_' side(1)],close2centr};
                            if ~ismember(lfp_act_chan,mychan)
                                continue
                            else
                                fprintf('Close 2 %s centr %s %s: %s\n',close2centr,patnm,side,lfp_act_chan);
                            end
                            
                        elseif ~strcmp('all',maxpower)
                            load(fullfile(rootdir,'maxFREQ_channels.mat'));
                            mychan =  maxFREQ_channels{[patnm '_' side(1)],maxpower};
                            if ~ismember(lfp_act_chan,mychan)
                                continue
                            else
                                fprintf('Channel with max %s power, %s %s: %s\n',maxpower,patnm,side,lfp_act_chan);
                            end
                            
                        end
                        
                        if isempty(lfpavg)
                            lfpavg{1} = EventAVGs.(lfp_act_chan).(basnm);
                        else
                            lfpavg{end+1} = EventAVGs.(lfp_act_chan).(basnm);
                        end
                    end
                    
                    if isempty(AP.(act_chan).AvgPow.(event).(subevent))
                        AP.(act_chan).AvgPow.(event).(subevent){1} = mean(cat(3,lfpavg{:}),3);
                    else
                        AP.(act_chan).AvgPow.(event).(subevent){end+1} = mean(cat(3,lfpavg{:}),3);
                    end
                end
                
            end
            
        end
        
    end
end
end




%--------------------------------------------------------------------------
function ssrt_topoplot(sess2analyse,varargin)

% Creates topoplots of postop EEG data data.
% Uses results average time freq maps (see save_avgepoch_pow).
% Calculates average across a specific frequency band and time windows of
% 250 ms. The figure includes data of [-1 1] sec relative of event (= 8
% topoplots for one group condition) and two group conditions (= 16
% topoplots in total).

% varargin = {};

prs = inputParser;
addRequired(prs,'sess2analyse',@isstruct);
% addParameter(prs,'grcond_names',{'left_stimoff','left_stimon','right_stimoff','right_stimon'},@iscell);
addParameter(prs,'topo_freq_nms',{'delta','theta','beta'},@iscell);
addParameter(prs,'topo_freqs',[0.5 4; 4, 7; 13,30],@ismatrix);
addParameter(prs,'EventTypes',{'StimulusOn','StopSignal'},@iscell);
addParameter(prs,'SubEventTypes', {'FailedStopTrial','SuccesfulStopTrial';...
    'FailedStopTrial','SuccesfulStopTrial'},@iscell);
addParameter(prs,'rectype','EEG',@ischar);
addParameter(prs,'csd',false,@islogical);
addParameter(prs,'epoch_win',[-2 2],@isvector);
addParameter(prs,'plot_win',[-0.5 0.5],@isvector);
addParameter(prs,'topobin',250,@isnumeric);
addParameter(prs,'sr',250,@isnumeric);
addParameter(prs,'cLim',[0.5 1.1],@(x) isvector(x)|isnumeric(x));
addParameter(prs,'preproc','epochs',@ischar);
addParameter(prs,'measure','Pow',@ischar);
addParameter(prs,'individual',true,@islogical);
addParameter(prs,'patavg',true,@islogical);
addParameter(prs,'grcond_names',{'left_stimoff'},@iscell);
addParameter(prs,'stat',true,@islogical);
addParameter(prs,'baseline_win',[-1 -.5],@(x) isvector(x)|isnumeric(x));
parse(prs,sess2analyse,varargin{:});
p = prs.Results;

global rootdir filesdir figdir_pd


%%


tfdir0 =fullfile(figdir_pd,'postop_EEG','TF');
topodir = fullfile(figdir_pd,'Postop_EEG','Topoplots');



topodir_rt= fullfile([topodir '_CSD' char(string(p.csd))]);
tfdir =[tfdir0 '_CSD' char(string(p.csd))];
res_nm = ['avgepoch_CSD' char(string(p.csd))];


chanlocs = [];
load(fullfile(rootdir,'postop_EEG_chanlocs.mat'));
chanlabels = {chanlocs.labels};




% Time windows

ds_sr = 50;
epoch_times =  p.epoch_win(1):1/ds_sr:p.epoch_win(2);

plotlims = dsearchn(epoch_times',p.plot_win');

if ~isempty(p.baseline_win)
    baslims = dsearchn(epoch_times',p.baseline_win');
    
    basnm = ['BASLims' num2str(baslims(1)) '_' num2str(baslims(2))]; basnm2 = basnm;
    plotlims = plotlims-baslims(1)+1;
    plotinx = plotlims(1):plotlims(2);
else
    basnm = 'TRN_AVG';
    basnm2 = 'TrialNorm';
    
    plotinx = plotlims(1):plotlims(2);
end

if isempty(p.SubEventTypes); subevs = 3; else; subevs = 1:3; end;

% Load freq components
f = [];
load(fullfile(rootdir,'freq_components.mat'))
%% Topoplots for individual patients (averaged across epochs)
if p.individual
    
    
    for snr = 1:length(sess2analyse)
        curr_resdir = sess2analyse(snr).folder;
        patnm = sess2analyse(snr).patient;
        side  = sess2analyse(snr).side;
        condition  = sess2analyse(snr).tag;
        
        resdir = fullfile(curr_resdir, ['EventAVGs_CSD' char(string(p.csd))]);
        
        
        for ei = 1:length(p.EventTypes)
            event = p.EventTypes{ei};
            for sei = subevs
                if sei<3
                    subevent = p.SubEventTypes{ei,sei};
                else
                    subevent = event;
                end
                
                
                evavg_nm = [ event '_' subevent '_AVGs.mat'];
                evSTAT_nm = [ event '_' subevent '_STATs.mat'];
                
                try
                    load(fullfile(resdir,evavg_nm));
                catch
                    fprintf('No EventAVG for %s\n',curr_resdir)
                    continue
                end
                
                if p.stat
                    load(fullfile(resdir,evSTAT_nm));
                end
                
                % Load patient data
                chanlabs = fieldnames(EventAVGs);
                alldat_epoch = cell(1,length(chanlabs));
                for ch = 1:length(chanlabs)
                    act_chan = chanlabs{ch};
                    alldat_epoch{ch} = EventAVGs.(act_chan).(basnm);
                end
                
                alldat_epoch2 = cat(3,alldat_epoch{:});
                
                
                
                % Make plots for each frequency band
                for fri = 1:length(p.topo_freq_nms)
                    
                    freqs = p.topo_freqs(fri,:);
                    
                    f_ind = intersect(find(f>=freqs(1)),find(f<=freqs(2)));
                    
                    
                    
                    % Select data range to plot
                    alldat = alldat_epoch2(f_ind,plotinx,:);
                    
                    % Stat
                    if p.stat
                        mask_ersp = cell(1,length(chanlabs)); trn = nan(1,length(chanlabs));
                        for ch = 1:length(chanlabs)
                            act_chan = chanlabs{ch};
                            mask_ersp{ch} = EventSTAT.(act_chan).(basnm2).mask_ersp;
                            trn(ch) = EventSTAT.(act_chan).(basnm2).trialnr;
                        end
                        mask_ersp_all = cat(3,mask_ersp{:});
                        mask_ersp_datwin = mask_ersp_all(f_ind,plotinx,:);
                        trnr = max(trn);
                        pattit = [patnm ' ' side ' ' condition '(n=' num2str(trnr) 'trials) '];
                    else
                        mask_topo = [];
                        mask_ersp_datwin = [];
                        pattit = [patnm ' ' side ' ' condition];
                    end
                    suptit = {[num2str(freqs(1)) '-' num2str(freqs(2)) 'Hz'],pattit,event, subevent};
                    
                    % Result directory to save figures
                    topodir = fullfile(topodir_rt,[event '_' subevent],...
                        ['WIN' num2str(p.plot_win(1)) '_' num2str(p.plot_win(2)) '_BIN' num2str(p.topobin) '_BAS' num2str(p.baseline_win)],p.topo_freq_nms{fri});
                    if ~isdir(topodir); mkdir(topodir); end;
                    
                    % TOPOPLOT
                    topoplot_fig(alldat,chanlocs,ds_sr,p.plot_win,p.topobin,...
                        p.cLim,topodir,mask_ersp_datwin,suptit,[patnm '_' side '_' condition])
                    
                end
            end
        end
    end
end

%%
if p.patavg
    
    
    for gri = 1:length(p.grcond_names)
        grnm = p.grcond_names{gri};
        grnm_tit = grnm; grnm_tit(strfind(grnm,'_')) = ' ';
        
        
        for ei = 1:length(p.EventTypes)
            event = p.EventTypes{ei};
            for sei = subevs
                if sei<3
                    subevent = p.SubEventTypes{ei,sei};
                else
                    subevent = event;
                end
                
                
                % Load group average
                alldat_epoch = [];
                for ch = 1:length(chanlabels)
                    act_chan = chanlabels{ch};
                    resdir_chan = fullfile(tfdir,act_chan);
                    
                    
                    
                    try
                        
                        load(fullfile(resdir_chan,[event '_' subevent '_' res_nm '_pow_' grnm  num2str(p.baseline_win) '.mat']));
                    catch
                        fprintf('No PatAVG for %s\n',resdir_chan);
                        continue
                    end
                    
                    if isempty(alldat_epoch)
                        alldat_epoch = nan(size(avgepoch_pow,1),size(avgepoch_pow,2),length(chanlabels));
                    end
                    
                    alldat_epoch(:,:,ch) = avgepoch_pow; avgepoch_pow = [];
                    
                end
                
                
                % Make plots for each frequency band
                for fri = 1:length(p.topo_freq_nms)
                    
                    freqs = p.topo_freqs(fri,:);
                    f_ind = intersect(find(f>=freqs(1)),find(f<=freqs(2)));
                    
                    
                    
                    fprintf('label: %s, side: %s, condition: %s, freq: %d-%d Hz (%s)\n',grnm, sess2analyse(1).side,  sess2analyse(1).tag, freqs,p.topo_freq_nms{fri})
                    
                    % Select data range to plot
                    alldat = alldat_epoch(f_ind,plotinx,:);
                    
                    
                    
                    % Stat
                    if p.stat
                        
                        mask_ersp = cell(1,length(chanlabels)); patnr = nan(1,length(chanlabels));
                        for ch = 1:length(chanlabels)
                            
                            act_chan = chanlabels{ch};
                            resdir_chan = fullfile(tfdir,act_chan);
                            
                            statfnm = fullfile(resdir_chan,[event '_' subevent '_' res_nm '_STAT_' grnm num2str(p.baseline_win) '.mat']);
                            load(statfnm);
                            
                            
                            mask_ersp{ch} = avgepoch_STAT.mask_ersp;
                            patnr(ch) = avgepoch_STAT.patientnr;
                            
                            avgepoch_STAT = [];
                        end
                        mask_ersp_all = cat(3,mask_ersp{:});
                        mask_ersp_datwin = mask_ersp_all(f_ind,plotinx,:);
                        
                        patientnr = max(patnr);
                        
                        pattit = ['Patient AVG:' grnm_tit '(n=' num2str(patientnr) ')'];
                    else
                        mask_topo = [];
                        patientnr = [];
                        mask_ersp_datwin = [];
                        pattit = ['Patient AVG:' grnm_tit];
                    end
                    
                    try
                        
                        for ch = 1:length(chanlabels)
                            
                            act_chan = chanlabels{ch};
                            resdir_chan = fullfile(tfdir,act_chan);
                            statfnm = fullfile(resdir_chan,[event '_' subevent '_' res_nm '_STAT_' grnm num2str(p.baseline_win) '.mat']);
                            load(statfnm);
                            
                            patnr(ch) = avgepoch_STAT.patientnr;
                        end
                        patientnr = max(patnr);
                        pattit = ['Patient AVG:' grnm_tit '(n=' num2str(patientnr) ')'];
                    end
                    
                    
                    
                    suptit = {[num2str(freqs(1)) '-' num2str(freqs(2)) 'Hz'],pattit,event, subevent};
                    
                    % Result directory to save figures
                    topodir = fullfile(topodir_rt,[event '_' subevent],...
                        ['WIN' num2str(p.plot_win(1)) '_' num2str(p.plot_win(2)) '_BIN' num2str(p.topobin) '_BAS' num2str(p.baseline_win)],p.topo_freq_nms{fri} );
                    if ~isdir(topodir); mkdir(topodir); end;
                    
                    % TOPOPLOT
                    
                    
                    topoplot_fig(alldat,chanlocs,ds_sr,p.plot_win,p.topobin,...
                        p.cLim,topodir,mask_ersp_datwin,suptit,grnm)
                end
            end
        end
    end
end
end
