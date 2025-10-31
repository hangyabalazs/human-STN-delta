function PD_eeg_stats(sess2analyse,EventTypes,SubEventTypes,groups2comp_labels,varargin)
%PD_EEG_STATS   Comparison of time-frequency maps
%   PD_eeg_stats(sess2analyse,EventTypes,SubEventTypes,groups2comp_labels,...)
%         performs intra and inter subject comparison of specified event partitions/ patient groups/ task
%         conditions using permutation test with or without correction for
%         multiple comparison (using std_data.m by EEGLAB or statcondfieldtrip.m by Fieldtrip toolboxes).
%         Plots group/ condition ERSP averages along with ERSP difference
%         maps. Topoplot averages and - difference maps also plotted for
%         postoperative EEG data.
%         Handles one or two  pairs of conditions/ groups.
%

% Required inputs:
%     SESS2ANALYSE          struct containing all necessary information (name of patient, side
%                           of experiment, tag of condition, session folder path) of
%                           session data that need to be analysed (see getdata2analyse)
%
%     EVENTTYPES            1xN cell array of event labels, ex: {'StimulusOn','StopSignal','KeyPress1','Feedback'};
%
%     SUBEVENTTYPES         Nx2 cell array of partition ("subevent") labels, each row
%                           corresponds to an event label, columns to partitions
%                           {'CueStim','StopStim';'FailedStopTrial','SuccesfulStopTrial';'CueResponse','StopResponse';'Correct','Error';};
%
%     GROUPS2COMP_LABELS    labels of groups/ conditions to compare
%           {'conditions'} | {'partitions'} | {'patgroups'} | {'conditions','partitions'} | {'conditions','patgroups'}
%
%
% Optional inputs (name-value pairs with default value):
%   'patgroup_nm'   cell array of patient group labels
%               {}, includes all patients
%               {'tremor-dominant','akinetic-rigid','mixed'}, draws separate plots for each clinical group
%               {'RTdecrease','RTincrease'}, - draws separate plots for each patient group based on preop-postop RT change
%               (default value: {})
%
%   'fr_names'      cell array, list of frequency band labels to compare
%                   (default value: {'delta'})
%
%   'freq_bands'    matrix, frequency band limits, each row corresponds to
%                   a label in 'fr_names' (def. value: [1 4]);
%
%   'topo_wins'     matrix, first and last value respresents the time limits of ERSPs, 
%                   if multiple rows, each row  corresponds to the time limits of one topoplot
%                   (power data averaged within this time bin) (def. value: [-1 -0.5; -0.5 0; 0 0.5; 0.5 1])
%
%   'avg_clim'      color axis limit for power averages (def. value: [-1 1])
%
%   'diff_clim'     color axis limit for difference maps def. value: [-2 2])
%
%   'side'          char. array, tested side def. value: 'left')
%
%   'indivfig'      true | false, if true generates a plot for each channel
%                   (def. value: true)
%
%   'isfig'         true | false, if true figures are generated (def. value: true)
%
%   'what2run'      char. array
%           'ersp'          generates ERSP maps
%           'timeseries'    generates time series plot (averaged across frequency
%                           components of ERSP
%           'both'          ERSP + timeseries (def. value: {'ersp','timeseries'})
%
%   'csd'           true | false, if true, CSD transformed EEG data is used
%                   (relevant only for postoperative EEG data) (default value: true)
%
%   'bipol'         true | false, if true, F4-F3 bipolar derivation of EEG data is used
%                   (relevant only for postoperative EEG data) (default value: false)
%
%   'chanmean'      if true, channel averaged LFP data is used, if false
%                   LFP data is plotted channel-by-channel
%                   (relevant only for intraop LFP data)
%                   1 | 0 (default value: 1)
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
%   'baseline_win'  1x2 vector, time window relative to event timestamp in
%   sec, for baseline correction (def. value: [-1 -.5]);
%
%   'baseline_type' char. array, type of baseline correction to apply
%           'common' | 'indiv' (def. value: 'common')
%
%   'alphas'        numeric or vector, list of alpha threshold to consider when performing statistics
%                   ERSP masks (1 where stat. significant, 0 otherwise) are created using the the specified alpha threshold
%                   if NaN, no statistical test is performed
%                   (default value: [0.05])
%
%   'stat_type'     statistical algorighm to use
%           'eeglab' | 'fieldtrip' (def. value: 'fieldtrip')
%
%   'mcorrect'     method to correct for multiple comparisons
%           'fdr' | 'cluster' | 'none' (def. value: 'cluster')
%
%   'bar_win'      1x2 vector, time limits for boxplots to compare (average
%                   across time and freq) (def. value: [0 0.5])
%
%   'bar_freq'     1x2 vector, freq limits for boxplots to compare (def. value: [1 2])
%
%   'intrastat'    true | false, INTRAsubject comparison is performed (def. value: true)
%
%   'interstat'    true | false, , INTERsubject comparison is performed (def. value: true)
%
% See also: TIME_FREQ_PATIENTS, STD_STAT, STATCONDFIELDTRIP

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu



prs = inputParser;
addRequired(prs,'sess2analyse',@isstruct);
addRequired(prs,'EventTypes',@iscell);
addRequired(prs,'SubEventTypes',@iscell);
addRequired(prs,'groups2comp_labels',@iscell);
addParameter(prs,'patgroup_nm',{},@iscell);

addParameter(prs,'fr_names',{'wide_band'},@iscell);
addParameter(prs,'freq_bands',[1 80],@ismatrix);

addParameter(prs,'topo_wins',[-1 -0.5; -0.5 0; 0 0.5; 0.5 1],@ismatrix);
addParameter(prs,'avg_clim',[-1 1],@isvector);
addParameter(prs,'diff_clim',[-2 2],@isvector);
addParameter(prs,'side','left',@ischar);
addParameter(prs,'indivfig',true,@islogical);
addParameter(prs,'isfig',true,@islogical);
addParameter(prs,'what2run',{'ersp','timeseries'},@(x) ischar(x)|iscell(x)); % 'ersp' | 'timeseries' | 'bar' | {'ersp','timeseries'}
addParameter(prs,'csd',false,@islogical);
addParameter(prs,'bipol',false,@islogical);
addParameter(prs,'measure','Pow',@ischar);

addParameter(prs,'chanmean',1,@isnumeric);
addParameter(prs,'subregion','all',@(x) ischar(x)||iscell(x));
addParameter(prs,'close2centr','all',@ischar);
addParameter(prs,'maxpower','all',@ischar);

addParameter(prs,'baseline_win',[-1 -.5],@(x) isnumeric(x)||isvector(x));
addParameter(prs,'baseline_type','common',@ischar);

addParameter(prs,'alphas',[0.05],@(x) isnumeric(x)||isvector(x));
addParameter(prs,'stat_type','fieldtrip',@ischar);
addParameter(prs,'mcorrect','cluster',@ischar);
addParameter(prs,'bar_win',[0 0.5],@isvector);
addParameter(prs,'bar_freq',[1 2],@isvector);


addParameter(prs,'intrastat',true,@islogical);
addParameter(prs,'interstat',true,@islogical);

parse(prs,sess2analyse,EventTypes,SubEventTypes,groups2comp_labels,varargin{:});
p = prs.Results;



compnr = length(groups2comp_labels);
rectime = sess2analyse(1).rectime;
rectype = sess2analyse(1).rectype;
tag = unique({sess2analyse.tag});
if length(tag)==1; cond = tag{1}; 
else; cond = 'bothcond'; 
end;

% %%% Within subject comparison %%%
% if p.intrastat
%     statwin = [p.topo_wins(1,1) p.topo_wins(end,end)];
%     
%     for alfi = p.alphas
%         within_pat_stat_loop(sess2analyse,EventTypes,SubEventTypes,p.freq_bands,...
%             statwin,p.subregion,alfi,p.mcorrect,p.avg_clim,p.diff_clim,p.baseline_win,p.csd);
%     end
% end

%%% Inter subject comparison %%%
if p.interstat
    % Select sessions according to side
    if ~contains(p.side,'both')
        sidinx = ismember({sess2analyse.side},p.side);
        sess2analyse2 = sess2analyse(sidinx);
    else
        sess2analyse2 = sess2analyse;
    end
    
    
    
    % Sort sessions according to clinical group
    sess2analyse3 = {};
    if ~isempty(p.patgroup_nm)
        patgroups = clinical_groups(p.patgroup_nm,rectime,cond);
        if isempty(patgroups)
            fprintf('No patient in this group\n'); return;
        end
        for k = 1:length(patgroups)
            patinx = ismember({sess2analyse2.patient},patgroups{k});
            sess2analyse3{k} = sess2analyse2(patinx);
        end
        sessgroup =  [p.side '_RACEfalse_' cat(2,p.patgroup_nm{:})];
    else
        sess2analyse3{1} = sess2analyse2;
        sessgroup = [p.side '_RACEfalse'];
    end
    
    if contains(p.side,'both')||length(tag)==1
        sessgroup = [sessgroup '_' cond];
    end
    %% Frequency
    for fii = 1:length(p.fr_names)
        frx = p.freq_bands(fii,:);
        frnm = p.fr_names{fii};
        
        
        %% Event
        for ei = 1:length(EventTypes)
            
            event2comp = EventTypes{ei};
            
            
            %% Groups to compare
            if ismember('partitions',groups2comp_labels) && compnr==1
                groups2comp = SubEventTypes(ei,:);
                
            elseif ismember('conditions',groups2comp_labels) && compnr==1
                groups2comp = {'stimoff','stimon'};
                
            elseif ismember('patgroups',groups2comp_labels) && compnr==1
                groups2comp = {p.patgroup_nm{:}};
                
            elseif ismember('partitions',groups2comp_labels) && ismember('conditions',groups2comp_labels)
                groups2comp = {'stimoff','stimon'; SubEventTypes{ei,:}};
                
            elseif ismember('patgroups',groups2comp_labels) && ismember('conditions',groups2comp_labels)
                groups2comp = {'stimoff','stimon'; p.patgroup_nm{:}};
                
            elseif ismember('sides',groups2comp_labels) && compnr==1
                groups2comp = {'left','right'};
                
                
                
            end
            
            
            %% Alpha level
            for alfi = p.alphas
                
                if contains('ersp',p.what2run)
                    %% STAT
                    
                    stat_twoway(sess2analyse3,groups2comp_labels,groups2comp,event2comp,...
                        'freq_bands',frx,'fr_name',frnm,'stat',p.stat_type,'alpha',alfi, ...
                        'topo_wins',p.topo_wins,'indivfig',p.indivfig,'isfig',p.isfig,...
                        'subregion',p.subregion,'close2centr',p.close2centr,'maxpower',p.maxpower,...
                        'dominantfreq',true,'sessgroup',sessgroup,'csd',p.csd,'bipol',p.bipol,...
                        'measure',p.measure,'avg_clim',p.avg_clim,'diff_clim',p.diff_clim,...
                        'rectime',rectime,'rectype',rectype,'side',p.side,'mcorrect',p.mcorrect,...
                        'baseline_win',p.baseline_win,'baseline_type', p.baseline_type,'chanmean',p.chanmean);
                    
                end
                
                
                
                
                if contains('timeseries',p.what2run)
                    
                    % Power time series from using statmat saved by stat_twoway
                    pow_timeseries_statmat(groups2comp_labels,groups2comp,event2comp,...
                        'fr_name',frnm,'stat',p.stat_type,'alpha',alfi,'csd',p.csd,...
                        'topo_wins',p.topo_wins,'sessgroup',sessgroup,...
                        'rectime',rectime,'rectype',rectype,'side',p.side,...
                        'avg_clim',p.avg_clim,'errorshade','SE','baseline_type',p.baseline_type);
                end
                
                
                if contains('bar',p.what2run)
                    
                    bargraphs_statmat(p.bar_win,p.bar_freq, groups2comp_labels,groups2comp,event2comp,...
                        'fr_name',frnm,'freq_bands',frx,'stat',p.stat_type,'alpha',alfi, ...
                        'topo_wins',p.topo_wins,'sessgroup',sessgroup,'csd',p.csd,...
                        'rectime',rectime,'rectype',rectype,'side',p.side)
                end
                
                
            end
            
        end
    end
end
end


%--------------------------------------------------------------------------

function stat_twoway(sess2analyse,groups2comp_labels,groups2comp,event2comp,varargin)

dbstop if error

prs = inputParser;
addRequired(prs,'sess2analyse',@(x) isstruct(x)||iscell(x));

addRequired(prs,'groups2comp_labels',@iscell);
addRequired(prs,'groups2comp',@iscell);
addRequired(prs,'event2comp',@ischar);
addParameter(prs,'sessgroup','left',@ischar);

addParameter(prs,'topo_wins',[-0.5 0; 0 0.5],@ismatrix);
addParameter(prs,'stat','eeglab',@ischar);
addParameter(prs,'alpha',0.01,@isnumeric);
addParameter(prs,'mcorrect','fdr',@ischar);
addParameter(prs,'baseline_win',[],@(x) isvector(x)|isnumeric(x))
addParameter(prs,'baseline_type','common',@ischar)

addParameter(prs,'freq_bands',[1 4],@isvector);
addParameter(prs,'fr_name','delta',@isvector);
addParameter(prs,'dominantfreq',false,@islogical);

addParameter(prs,'indivfig',true,@islogical);
addParameter(prs,'isfig',true,@islogical);
addParameter(prs,'avg_clim',[0.5 1.2],@isvector);
addParameter(prs,'diff_clim',[-2 2],@isvector);

addParameter(prs,'chanmean',1,@isnumeric);
addParameter(prs,'subregion','all',@ischar);
addParameter(prs,'close2centr','all',@ischar);
addParameter(prs,'maxpower','all',@ischar);

addParameter(prs,'preproc','epochs',@ischar);
addParameter(prs,'csd',false,@islogical);
addParameter(prs,'bipol',false,@islogical);
addParameter(prs,'measure','Pow',@ischar);
addParameter(prs,'rectime','postop',@ischar);
addParameter(prs,'rectype','EEG',@ischar);
addParameter(prs,'side','left',@ischar);
parse(prs,sess2analyse,groups2comp_labels,groups2comp,event2comp,varargin{:})
p = prs.Results;



global rootdir figdir_pd

%% Parameters & directories



if strcmp(p.rectime,'postop')
    chanlocs = [];
    if ~p.bipol
        load(fullfile(rootdir,'postop_EEG_chanlocs.mat'));
        chanlabels = {chanlocs.labels};
    else
        chanlabels = {'F4_F3'};
    end
elseif strcmp(p.rectime,'intraop') && strcmp(p.rectype,'EEG')
    chanlabels = {'F4'};
elseif strcmp(p.rectime,'intraop') && strcmp(p.rectype,'LFP')
    chanlabels = {'chanmean'};
end
chnr = length(chanlabels);


%% Label for groups to compare (for saving)
compnr = length(p.groups2comp_labels);

if compnr==1
    
    gr2c = [p.groups2comp{1} '_vs_' p.groups2comp{2}];
else
    
    gr2c = [p.groups2comp{1,1} '_vs_' p.groups2comp{1,2} '_' p.groups2comp{2,1} '_' p.groups2comp{2,2}];
end

if ismember('partitions',p.groups2comp_labels)
    partnr = 2;
    subevent2comp = p.groups2comp(ismember(p.groups2comp_labels,'partitions'),:);
else
    partnr = 1;
    subevent2comp = {p.event2comp};
end



%% Time


% Plot windows
wnnm = num2str(p.topo_wins);
wnnm = wnnm(:)';
stat_win = [p.topo_wins(1,1) p.topo_wins(end,end)];


epoch_win = [-2 2];
new_sr = 50;
times = epoch_win(1):1/new_sr:epoch_win(2);
new_times = stat_win(1):1/new_sr:stat_win(2);
newtinx = dsearchn(times',new_times');

% Baseline window
if ~isempty(p.baseline_win)
    baslims = dsearchn(times',p.baseline_win');
    
    
    if strcmp(p.baseline_type,'indiv')
        newtinx = newtinx-baslims(1)+1;
    end
    
    
end

%% Result directory

resdir_statfigs = fullfile(figdir_pd,[p.rectime '_' p.rectype],'stats',p.sessgroup , ...
    p.fr_name,p.event2comp,gr2c,['CSD' char(string(p.csd)) '_' p.stat '_Win' wnnm]);
if ~isdir(resdir_statfigs);mkdir(resdir_statfigs); end;


%% Frequency components derived from wavelet transform
load(fullfile(rootdir,'freq_components.mat'));
f_ind = intersect(find(f>=p.freq_bands(1)),find(f<=p.freq_bands(2)));


%% Localization
if ~strcmp(p.subregion,'all') && strcmp(p.rectype,'LFP')
    load(fullfile(rootdir,'STN_loc','STN_loc.mat'));
    subreg_names = STN_loc.subreg_names;
    sinx = find(strcmp(subreg_names,p.subregion));
end





if exist(fullfile(resdir_statfigs,['statmatx_' p.side '_' num2str(p.alpha) '_' p.baseline_type '.mat']))~=2
    
    %% Patientgroups
    grnr = length(p.sess2analyse);
    
    grmat = cell(1,grnr);
    for gri = 1:grnr
        
        s2a = p.sess2analyse{gri};
        patnr = length(s2a);
        patmat = cell(patnr,chnr);
        
        %% Load data of all patients (from 1 patient group)
        for ip = 1:patnr
            
            
            patnm = s2a(ip).patient;
            curr_resdir = s2a(ip).folder;
            
            if ~p.bipol
                resdir = fullfile(curr_resdir, ['EventAVGs_CSD' char(string(p.csd))]);
            else
                resdir = fullfile(curr_resdir,'EventAVGs_bipol');
            end
            
%             Individual baseline (/ partition/ one patient)
            if strcmp(p.baseline_type,'indiv') && ~isempty(p.baseline_win)
                basnm = ['BASLims' num2str(baslims(1)) '_' num2str(baslims(2))];
            else
                basnm = ['TRN_AVG'];
            end
            
            
            
            for ich = 1:chnr
                
                if ~strcmp(p.rectype,'LFP') || (strcmp(p.rectype,'LFP') && p.chanmean==1)
                    act_chan = chanlabels{ich};
                    
                    
                    % Loop over partitions (1 - Failed Stop, 2 - Successful Stop)
                    for k = 1:partnr
                        
                        if ~strcmp(p.rectype,'LFP')
                            evavg_nm = [ p.event2comp '_' subevent2comp{k} '_AVGs.mat'];
                            
                        elseif (strcmp(p.rectype,'LFP') && p.chanmean==1)
                            evavg_nm = [ p.event2comp '_' subevent2comp{k} '_chan_AVGs.mat'];
                        end
                        
                        try
                            load(fullfile(resdir,evavg_nm));
                        catch
                            fprintf('No %s trial. %s \n',subevent2comp{k}, curr_resdir);
                            continue
                        end
                        %
                        try
                            patmat{ip,ich,k} = EventAVGs.(act_chan).(basnm);
                        catch
                            fprintf('No %s bascorr. %s \n',subevent2comp{k}, curr_resdir);
                        end
                        
                    end
                    
                    
                    
                    %                     end
                    
                    
                    
                    
                else
                    
                    origchanlab = arrayfun(@(x) ['Ch' num2str(x)], 1:5, 'UniformOutput',0);
                    
                    for k = 1:partnr
                        
                        evavg_nm = [ p.event2comp '_' subevent2comp{k} '_AVGs.mat'];
                        
                        try
                            load(fullfile(resdir,evavg_nm));
                        catch
                            fprintf('No %s trial.\n',subevent2comp{k});
                            continue
                        end
                        
                        powavg = {};
                        for cici = 1:length(origchanlab)
                            
                            ach = origchanlab{cici};
                            
                            
                            if ~isfield(EventAVGs,ach)
                                fprintf('No %s channel.\n',ach);
                                continue
                            end
                            
                            
                            if contains(lower(p.subregion),'no')
                                if any(~ismember(STN_loc.Patients.(patnm).(p.side).Channels.(ach).Subregion,-1));
                                    fprintf('Not %s: %s %s %s\n',p.subregion,patnm,p.side,ach)
                                    continue
                                end
                            elseif ~strcmp(p.subregion,'all')  && strcmp(p.rectype,'LFP')
                                if ~ismember(STN_loc.Patients.(patnm).(p.side).Channels.(ach).Subregion(sinx),[1 2]);
                                    fprintf('Not %s: %s %s %s\n',p.subregion,patnm,p.side,ach)
                                    continue
                                end
                            elseif ~strcmp('all',p.close2centr) && strcmp(p.rectype,'LFP')
                                load(fullfile(rootdir,'Channels_close2centroids.mat'));
                                mychan = centr_tab{[patnm '_' p.side(1)],p.close2centr};
                                if ~ismember(ach,mychan)
                                    continue
                                else
                                    fprintf('Close 2 %s centr %s %s: %s\n',p.close2centr,patnm,p.side,ach);
                                end
                                
                            elseif ~strcmp('all',p.maxpower) && strcmp(p.rectype,'LFP')
                                load(fullfile(rootdir,'maxFREQ_channels.mat'));
                                mychan =  maxFREQ_channels{[patnm '_' p.side(1)],p.maxpower};
                                if ~ismember(ach,mychan)
                                    continue
                                else
                                    fprintf('Channel with max %s power, %s %s: %s\n',p.maxpower,patnm,p.side,ach);
                                end
                                
                                
                            end
                            
                            if isempty(powavg)
                                powavg{1} = EventAVGs.(ach).(basnm);
                            else
                                powavg{end+1} = EventAVGs.(ach).(basnm);
                            end
                            
                        end
                        
                        
                        patmat{ip,ich,k} = mean(cat(3,powavg{:}),3);
                    end
                    
                    
                end
                
                
                %  Apply common baseline  (/ one patient)
                if strcmp(p.baseline_type,'common') && partnr==2
                    parts2b = cat(3,patmat{ip,ich,:}); % concat partitions (averaged across trials; one channel of one patient)
                    if ~isempty(parts2b)
                        bas_avg = repmat(    nanmean(parts2b(:,baslims(1): baslims(2),:),[2 3])                       ,[1 size(parts2b,2) 1]); % average over baseline period & partitions
                        bas_sd = repmat(     std(   parts2b(:,baslims(1): baslims(2),:)     ,[],[2 3], 'omitnan')     ,[1 size(parts2b,2) 1] ); % standard dev
                        for kk = 1:partnr
                            if ~isempty(patmat{ip,ich,kk})
                                patmat{ip,ich,kk} = (patmat{ip,ich,kk}- bas_avg)./bas_sd; % baseline correction applied for each partition
                            end
                        end
                    end
                end
                
                
                
            end
            
        end
        
        
        
        
        
        
        grmat{gri} = patmat;
    end
    
    
    
    %% PatientAVG ERSP
    
    
    % TF1 = cat(3,patmat{:,1,1});
    % tf1 = mean(TF1(f_ind,newtinx,:),3);
    %
    %
    % TF2 = cat(3,patmat{:,1,2});
    % tf2 = mean(TF2(f_ind,newtinx,:),3);
    % TFTF = (tf1+tf2)/2;
    % [tf,pow2,im] = spectr_fig(TFTF,f(f_ind),new_sr,p.freq_bands(2),p.freq_bands(1),stat_win, [], [],[0.5 1.2])
    
    %% Create matrix for statistics (for one channel data)
    
    statmat = cell(2,chnr); % patientgroups x channels
    % if partitioned: each cell contains 1x2 cell array containing 2D matrices (freq x time);
    % if not partitioned: each cell contains one 2D matrix
    
    for ic = 1:chnr
        
        
        if  ismember('conditions',groups2comp_labels) && compnr==1 % compare stimoff- stimon conditions
            for gi = 1:2
                cnd = p.groups2comp{gi};
                
                pinx = find(strcmp({p.sess2analyse{1}.tag},cnd));
                pow = cat(3,patmat{pinx,ic,1});
                
                if strcmp(p.baseline_type,'indiv')
                    bas_avg = repmat( nanmean(pow(:,baslims(1): baslims(2),:),[2 3])   ,[1 size(pow,2) size(pow,3)]);
                    bas_sd = repmat( std(   pow(:,baslims(1): baslims(2),:)     ,[],[2 3], 'omitnan')   ,[1 size(pow,2) size(pow,3)]);
                    pow = (pow- bas_avg)./ bas_sd;
                end
                statmat{gi,ic} = pow(f_ind,newtinx,:);
                
            end
            
            
            
        elseif ismember('partitions',groups2comp_labels) && compnr==1 % compare partitions
            pow1 = cat(3,patmat{:,ic,1}); % freq x time x patients
            pow2 = cat(3,patmat{:,ic,2});
            
            statmat{1,ic} = pow1(f_ind,newtinx,:);
            statmat{2,ic} = pow2(f_ind,newtinx,:);
            
            
        elseif ismember('patgroups',groups2comp_labels) && compnr==1 % compare partitions
            
            
            patgroup = cell(1,2);
            for pg = 1:2
                patmat = grmat{pg};
                patgroup{pg} = cat(3,patmat{:,ic});
            end
            statmat{1,ic} = patgroup{1}(f_ind,newtinx,:);
            statmat{2,ic} = patgroup{2}(f_ind,newtinx,:);
            
            
            
        elseif ismember('partitions',groups2comp_labels) && ismember('conditions',groups2comp_labels)
            for gi = 1:2 % 1- stimoff, 2- stimon
                cnd = p.groups2comp{1,gi};
                
                pinx{gi} = find(strcmp({p.sess2analyse{1}.tag},cnd));
                part1 = cat(3,patmat{pinx{gi},ic,1}); % Failed
                part2 = cat(3,patmat{pinx{gi},ic,2}); % Successful
                statmat{gi,ic} = {part1(f_ind,newtinx,:), part2(f_ind,newtinx,:)};
                cnd_nr(gi) = size(part1,3);
                
            end
            
            
        elseif  ismember('patgroups',groups2comp_labels) && ismember('conditions',groups2comp_labels)
            
            for gi = 1:2 % 1- stimoff, 2- stimon
                cnd = p.groups2comp{1,gi};
                
                patgroup = cell(1,2);
                for pg = 1:2 % pat groups
                    pinx{gi} = find(strcmp({p.sess2analyse{pg}.tag},cnd));
                    patmat = grmat{pg};
                    patgroup{pg} = cat(3,patmat{pinx{gi},ic});
                end
                statmat{gi,ic} = {patgroup{1}(f_ind,newtinx,:), patgroup{2}(f_ind,newtinx,:)};
            end
            
            
        elseif  ismember('sides',groups2comp_labels) && compnr==1 % compare stimoff- stimon conditions
                for gi = 1:2
                    sd = p.groups2comp{gi};
                    
                    pinx = find(strcmp({p.sess2analyse{1}.side},sd));
                    pow = cat(3,patmat{pinx,ic,1});
                    
                    if strcmp(p.baseline_type,'indiv')
                        bas_avg = repmat( nanmean(pow(:,baslims(1): baslims(2),:),[2 3])   ,[1 size(pow,2) size(pow,3)]);
                        bas_sd = repmat( std(   pow(:,baslims(1): baslims(2),:)     ,[],[2 3], 'omitnan')   ,[1 size(pow,2) size(pow,3)]);
                        pow = (pow- bas_avg)./ bas_sd;
                    end
                    statmat{gi,ic} = pow(f_ind,newtinx,:);
                    
                end
                
                
                
            
            
        end
        
        
        if compnr==1 && partnr==1
            if strcmp(p.baseline_type,'indiv')
                
                m1 = mean(statmat{1,ic},3);
                m2 =  mean(statmat{2,ic},3);
                parts2b = cat(3,m1,m2);
                baslims2 = baslims-newtinx(1)+1;
                bas_avg = nanmean(parts2b(:,baslims2(1): baslims2(2),:),[2 3]);
                bas_sd = std(   parts2b(:,baslims2(1): baslims2(2),:)     ,[],[2 3], 'omitnan');
                for gi = 1:2
                    statmat{gi,ic} =  (statmat{gi,ic}-repmat(bas_avg,[1 size(parts2b,2) size(statmat{gi,ic},3)]) )./  repmat(bas_sd,[1 size(parts2b,2) size(statmat{gi,ic},3)]);
                end
            elseif   strcmp(p.baseline_type,'common')
                parts2b = cat(3,statmat{:,ic});
                baslims2 = baslims-newtinx(1)+1;
                bas_avg = nanmean(parts2b(:,baslims2(1): baslims2(2),:),[2 3]);
                bas_sd = std(   parts2b(:,baslims2(1): baslims2(2),:)     ,[],[2 3], 'omitnan');
                for gi = 1:2
                    statmat{gi,ic} =  (statmat{gi,ic}-repmat(bas_avg,[1 size(parts2b,2) size(statmat{gi,ic},3)]) )./  repmat(bas_sd,[1 size(parts2b,2) size(statmat{gi,ic},3)]);
                end
            end
        end
        
        
    end
    
    
    
    
    
    %% STAT
    
    
    % conditions: SuccessfulStop vs Failed Stop
    % groups: stimoff vs stimon
    
    
    for ichh = 1:chnr
        
        act_chan = chanlabels{ichh};
        
        
        if compnr==1
            
            %% One-way stat
            if strcmp(p.stat,'eeglab')
                [pcond{1,ichh}, ~, ~, statscond{1,ichh}, ~, ~] = std_stat(statmat(:,ichh),'condstats','on','groupstats','off','mode','eeglab',...
                    'method','perm','mcorrect',p.mcorrect,'alpha',p.alpha);
            elseif strcmp(p.stat,'fieldtrip')
                
                [F, df, pval]  = statcondfieldtrip(statmat(:,ichh)', 'method','montecarlo','naccu', 1000,...
                    'mcorrect',p.mcorrect,'alpha',p.alpha);
                pcond{1,ichh} = {pval<p.alpha};
                statscond{1,ichh} = {F};
            end
            
            if p.isfig
                fig = figure;
                set(fig,'Visible','off')
                
                for si = 1:2
                    subplot(3,1,si);
                    tf = mean(statmat{si,ichh},3);
                    [~,pow2,im] = spectr_fig(tf,f(f_ind),new_sr,p.freq_bands(2),p.freq_bands(1),stat_win, [], 0,p.avg_clim);
                    title([p.groups2comp{si} ', n='  num2str( size(statmat{si,ichh},3) ) ]);
                    xL = xlim; xticks(round(linspace(xL(1),xL(end)-0.5,5))); xticklabels(arrayfun(@num2str,new_times(xticks),'UniformOutput',0));
                    xlabel('Time (s)');
                    
                end
                
                subplot(3,1,3)
                tf3 = cell2mat(statscond{1,ichh});
                [~,pow2,im] = spectr_fig(tf3,f(f_ind),new_sr,p.freq_bands(2),p.freq_bands(1),stat_win, [], 0,p.diff_clim);
                hold on;
                contour(cell2mat(pcond{1,ichh}),'Color','white'); title('Difference (F-S)')
                
                xL = xlim; xticks(round(linspace(xL(1),xL(end)-0.5,5))); xticklabels(arrayfun(@num2str,new_times(xticks),'UniformOutput',0));
                xlabel('Time (s)');
                
                suptitle({['Aligned to ' p.event2comp],act_chan,p.fr_name})
                set(fig,'Position',get(0,'Screensize'))
                
                if strcmp(p.rectype,'LFP')
                    fnm = fullfile(resdir_statfigs,[act_chan  '_' num2str(p.alpha) '_' p.subregion 'STN' '_CLOSE2' p.close2centr '_MAXPOW' p.maxpower '_' p.mcorrect '_BAS' num2str(p.baseline_win)]);
                else
                    fnm = fullfile(resdir_statfigs,[act_chan  '_' num2str(p.alpha) '_' p.mcorrect '_BAS' num2str(p.baseline_win)]);
                end
                fnm = [fnm '_' p.baseline_type];
                
                saveas(fig,[fnm '.jpg']);
                saveas(fig,[fnm '.fig']);
                saveas(fig,[fnm '.emf']); close(fig)
            end
            
            
        elseif compnr==2
            %% Two-way stat
            
            
            for st = 1:4
                
                switch st; case 1; STmat = cat(1,statmat{1,ichh})'; % stimoff Failed vs Succes
                    case 2; STmat = cat(1,statmat{2,ichh})'; % stimon Failed vs Succes
                    case 3; s1 =statmat{1,ichh}; s2 = statmat{2,ichh};
                        STmat = cat(1,s1(1), s2(1)); % stimoff Failed vs stimon Failed
                    case 4; s1 =statmat{1,ichh}; s2 = statmat{2,ichh};
                        STmat = cat(1,s1(2), s2(2)); end; % stimoff Succes vs stimon Succes
                    
                    
                    if strcmp(p.stat,'eeglab')
                        [pcond{st,ichh}, ~, ~, statscond{st,ichh}, ~, ~] = std_stat(STmat,'condstats','on','groupstats','off','mode','eeglab',...
                            'method','perm','mcorrect',p.mcorrect,'alpha',p.alpha);
                    elseif strcmp(p.stat,'fieldtrip')
                        [F, df, pval]  = statcondfieldtrip(STmat,'method','montecarlo','naccu',1000,...
                            'mcorrect',p.mcorrect,'alpha',p.alpha);
                        pcond{st,ichh} = {pval<p.alpha};
                        statscond{st,ichh} = {F};
                    end
            end
            
            if p.isfig
                if p.indivfig
                    
                    fig = figure;
                    set(fig,'Visible','off')
                    for cipi = 1:2
                        
                        for pipi = 1:2
                            spnr = (cipi-1)*3+pipi;
                            s = statmat{cipi,ichh};
                            
                            subplot(3,3,spnr);
                            tf =  mean(s{pipi},3);
                            [~,pow2,im] = spectr_fig(tf,f(f_ind),new_sr,p.freq_bands(2),p.freq_bands(1),stat_win, [], 0,p.avg_clim);
                            title({p.groups2comp{1,cipi},p.groups2comp{2,pipi}})
                            xL = xlim; xticks(round(linspace(xL(1),xL(end)-0.5,5))); xticklabels(arrayfun(@num2str,new_times(xticks),'UniformOutput',0));
                            xlabel('Time (s)');
                            
                            
                        end
                    end
                    
                    subplot(3,3,3)
                    tf =cell2mat(statscond{1,ichh});
                    [~,pow2,im] = spectr_fig(tf,f(f_ind),new_sr,p.freq_bands(2),p.freq_bands(1),stat_win, [], 0,p.diff_clim);
                    hold on;
                    contour(cell2mat(pcond{1,ichh}),'Color','white'); title('p values')
                    title([p.groups2comp{1,1} ': ' p.groups2comp{2,1} ' - ' p.groups2comp{2,2}])
                    xL = xlim; xticks(round(linspace(xL(1),xL(end)-0.5,5))); xticklabels(arrayfun(@num2str,new_times(xticks),'UniformOutput',0));
                    xlabel('Time (s)');
                    
                    
                    
                    subplot(3,3,6)
                    tf =cell2mat(statscond{2,ichh});
                    [~,pow2,im] = spectr_fig(tf,f(f_ind),new_sr,p.freq_bands(2),p.freq_bands(1),stat_win, [], 0,p.diff_clim);
                    hold on;
                    contour(cell2mat(pcond{2,ichh}),'Color','white'); title('p values')
                    title([p.groups2comp{1,2} ': ' p.groups2comp{2,1} ' - ' p.groups2comp{2,2}])
                    xL = xlim; xticks(round(linspace(xL(1),xL(end)-0.5,5))); xticklabels(arrayfun(@num2str,new_times(xticks),'UniformOutput',0));
                    xlabel('Time (s)');
                    
                    
                    subplot(3,3,7)
                    tf =cell2mat(statscond{3,ichh});
                    [~,pow2,im] = spectr_fig(tf,f(f_ind),new_sr,p.freq_bands(2),p.freq_bands(1),stat_win, [], 0,p.diff_clim);
                    hold on;
                    contour(cell2mat(pcond{3,ichh}),'Color','white'); title('p values')
                    title([p.groups2comp{2,1} ': ' p.groups2comp{1,1} ' - ' p.groups2comp{1,2}])
                    xL = xlim; xticks(round(linspace(xL(1),xL(end)-0.5,5))); xticklabels(arrayfun(@num2str,new_times(xticks),'UniformOutput',0));
                    xlabel('Time (s)');
                    
                    subplot(3,3,8)
                    tf =cell2mat(statscond{4,ichh});
                    [~,pow2,im] = spectr_fig(tf,f(f_ind),new_sr,p.freq_bands(2),p.freq_bands(1),stat_win, [], 0,p.diff_clim);
                    hold on;
                    contour(cell2mat(pcond{4,ichh}),'Color','white'); title('p values')
                    title([p.groups2comp{2,2} ': ' p.groups2comp{1,1} ' - ' p.groups2comp{1,2}])
                    xL = xlim; xticks(round(linspace(xL(1),xL(end)-0.5,5))); xticklabels(arrayfun(@num2str,new_times(xticks),'UniformOutput',0));
                    xlabel('Time (s)');
                    
                    
                    
                    set(0, 'DefaultFigureRenderer', 'painters');
                    suptitle({['Aligned to ' p.event2comp],act_chan,p.fr_name})
                    set(fig,'Position',get(0,'Screensize'))
                    
                    if strcmp(p.rectype,'LFP')
                        fnm = fullfile(resdir_statfigs,[p.side '_' act_chan  '_' num2str(p.alpha) '_' p.subregion 'STN_' p.mcorrect '_BAS' num2str(p.baseline_win)]);
                        
                    else
                        fnm = fullfile(resdir_statfigs,[p.side '_' act_chan  '_' num2str(p.alpha) '_' p.mcorrect '_BAS' num2str(p.baseline_win)]);
                    end
                    fnm = [fnm '_' p.baseline_type];
                    saveas(fig,[fnm '.jpg']);
                    saveas(fig,[fnm '.fig']); close(fig)
                end
                
                
            end
            
            
            
            
        end
        
    end
    
    if strcmp(p.rectype,'LFP')
        fnm = fullfile(resdir_statfigs,['statmatx_' p.side '_' num2str(p.alpha) '_' p.subregion 'STN_' p.baseline_type]);
    else
        fnm = fullfile(resdir_statfigs,['statmatx_' p.side '_' num2str(p.alpha) '_' p.baseline_type]);
    end
    fnm = [fnm '_' p.baseline_type];
    save([fnm '.mat'],'statmat','statscond','pcond','chanlabels');
else
    
    epoch_win = [-2 2];
    %     sr = EEG.srate;
    %     new_sr = round(200/diff(epoch_win));
    new_sr = 50;
    new_times = stat_win(1):1/new_sr:stat_win(2);
    
    
    load(fullfile(resdir_statfigs,['statmatx_' p.side '_' num2str(p.alpha) '_' p.baseline_type '.mat']))
    
end


%% Topoplots

if strcmp(p.rectime,'postop') && ~p.bipol
    
    winnr = size(p.topo_wins,1);
    topotinx = cell(winnr,1);
    for wi = 1:winnr
        topo_times = p.topo_wins(wi,1):1/new_sr:p.topo_wins(wi,2); topo_times = topo_times(1:end-1);
        topotinx{wi} = dsearchn(new_times',topo_times');
    end
    clear wi
    
    if p.isfig
        fig = figure;
        set(fig,'Visible','off')
        spcols = winnr*3+2;
        sprows = compnr+1;
        % Averages
        for ci = 1:2 % 1-stimoff 2-stimon
            
            spnr = (ci-1)*spcols+1;
            for pii = 1:compnr % 1-Failed 2-Successful
                
                if compnr==1
                    tf = statmat(ci,:);
                elseif compnr==2
                    s = statmat(ci,:);
                    tf = cellfun(@(x) x{pii},s,'UniformOutput',0);
                end
                
                
                for wiii = 1:winnr
                    subplot(sprows,spcols,spnr)
                    chvec1 = cell2mat(cellfun(@(x) mean(x(:,topotinx{wiii},:),'all'),tf,'UniformOutput',0));
                    
                    subjnr = size( tf{1} ,3);
                    
                    
                    topoplot(chvec1,chanlocs); caxis(p.avg_clim)
                    set_my_topo(gcf)
                    
                    
                    if wiii ==1
                        if compnr==1
                            title({[p.groups2comp{1,ci} ', n=' num2str(subjnr)],...
                                [num2str(p.topo_wins(wiii,1)) '  ' num2str(p.topo_wins(wiii,2)) ' s']});
                        elseif compnr==2
                            title({p.groups2comp{1,ci},[p.groups2comp{2,pii} ', n=' num2str(subjnr)],...
                                [num2str(p.topo_wins(wiii,1)) '  ' num2str(p.topo_wins(wiii,2)) ' s']});
                        end
                    else
                        
                        title({'','',[num2str(p.topo_wins(wiii,1)) '  ' num2str(p.topo_wins(wiii,2)) ' s']});
                    end
                    spnr = spnr+1;
                end
                
                spnr = spnr+1;
            end
        end
        
        
        % Difference maps
        statnr = size(pcond,1);
        for sti = 1:statnr
            switch sti; case 1; spnr = spcols-winnr+1; case 2; spnr = spcols*2-winnr+1;
                case 3; spnr =  spcols*2+1; case 4; spnr =  spcols*2+winnr+2; end;
                
                
                tf = cellfun(@(x) cell2mat(x),statscond(sti,:),'UniformOutput',0);
                pf = cellfun(@(x) cell2mat(x),pcond(sti,:),'UniformOutput',0);
                
                
                for wiii = 1:winnr
                    sp(spnr) = subplot(sprows,spcols,spnr);
                    chvec3 = cell2mat(cellfun(@(x) mean(x(:,topotinx{wiii}),'all'),tf,'UniformOutput',0));
                    
                    pvec3 = cell2mat(cellfun(@(x) any(find(x(:,topotinx{wiii}))),pf,'UniformOutput',0));
                    
                    topoplot(chvec3,chanlocs); caxis(p.diff_clim)
                    set_my_topo(gcf)
                    
                    if ~isempty(pvec3)
                        hold on; scatter(sp(spnr).Children(1).XData(pvec3),sp(spnr).Children(1).YData(pvec3),...
                            30,'white','filled')
                    end
                    
                    
                    title({'','',[num2str(p.topo_wins(wiii,1)) '  ' num2str(p.topo_wins(wiii,2)) ' s']});
                    
                    
                    spnr = spnr+1;
                end
                
        end
        subplot(sprows,spcols,sprows*spcols-winnr+1);
        topoplot(chvec1,chanlocs); caxis(p.avg_clim); colorbar;
        title('AVG topoplots');cla;
        
        
        subplot(sprows,spcols,sprows*spcols-winnr+2);
        topoplot(chvec3,chanlocs); caxis(p.diff_clim); colorbar;
        title('DIFF topoplots');cla;
        
        
        suptitle({['Aligned to ' p.event2comp], [num2str(p.freq_bands) ' Hz'],p.fr_name})
        
        set(fig,'Position',get(0,'Screensize'))
        set(0, 'DefaultFigureRenderer', 'painters');
        fnm = fullfile(resdir_statfigs,['Topo_' p.side '_' num2str(p.alpha) '_' p.mcorrect '_' p.baseline_type]);
        saveas(fig,[fnm '.jpg']);
        saveas(fig,[fnm '.fig']); close(fig);
    end
end
end


%--------------------------------------------------------------------------
function pow_timeseries_statmat(groups2comp_labels,groups2comp,event2comp,varargin)

prs = inputParser;

addRequired(prs,'groups2comp_labels',@iscell);
addRequired(prs,'groups2comp',@iscell);
addRequired(prs,'event2comp',@ischar);
addParameter(prs,'sessgroup','left',@ischar);

addParameter(prs,'topo_wins',[-0.5 0; 0 0.5],@ismatrix);
addParameter(prs,'stat','eeglab',@ischar);
addParameter(prs,'alpha',0.01,@isnumeric);
addParameter(prs,'mcorrect','fdr',@ischar);
addParameter(prs,'subregion','all',@ischar);

addParameter(prs,'fr_name','delta',@isvector);
addParameter(prs,'avg_clim',[-2 2],@isvector);

addParameter(prs,'csd',false,@islogical);

addParameter(prs,'rectime','postop',@ischar);
addParameter(prs,'rectype','EEG',@ischar);
addParameter(prs,'side','left',@ischar);

addParameter(prs,'errorshade','SE',@ischar);
addParameter(prs,'baseline_type','common',@ischar);
parse(prs,groups2comp_labels,groups2comp,event2comp,varargin{:})
p = prs.Results;



global figdir_pd

wnnm = num2str(p.topo_wins);
wnnm = wnnm(:)';
stat_win = [p.topo_wins(1,1) p.topo_wins(end,end)];
new_sr = 50;
new_times = stat_win(1):1/new_sr:stat_win(2);

compnr = length(p.groups2comp_labels);
if compnr==1
    
    gr2c = [p.groups2comp{1} '_vs_' p.groups2comp{2}];
else
    
    gr2c = [p.groups2comp{1,1} '_vs_' p.groups2comp{1,2} '_' p.groups2comp{2,1} '_' p.groups2comp{2,2}];
end


resdir_statfigs = fullfile(figdir_pd,[p.rectime '_' p.rectype],'stats',p.sessgroup,...
    p.fr_name,p.event2comp,gr2c,['CSD' char(string(p.csd)) '_' p.stat '_Win' wnnm]);


if strcmp(p.rectype,'LFP')
    load(fullfile(resdir_statfigs,['statmatx_' p.side '_' num2str(p.alpha) '_' p.subregion 'STN_' p.baseline_type '.mat']));
else
    load(fullfile(resdir_statfigs,['statmatx_' p.side '_' num2str(p.alpha) '_' p.baseline_type '.mat']));
end

channr = length(chanlabels);
for ichh = 1:channr
    fig = figure;
    set(fig,'Visible','off')
    
    colors = [1 0 0; 0 1 0; 0 0 0; 0 0 1];
    if compnr==1
        stval = 1;
    else
        stval = 1:4;
    end
    
    j = 1;
    for st = stval
        
        if compnr==1
            STmat = statmat(:,ichh);
            leg = {p.groups2comp{1,1}, p.groups2comp{1,2}};
            tit = '';
            cols = colors(1:2,:);
        else
            switch st; case 1; STmat = cat(1,statmat{1,ichh})'; % stimoff Failed vs Succes
                cols = colors(1:2,:);
                tit = p.groups2comp{1,1}; leg = {p.groups2comp{2,1}, p.groups2comp{2,2}};
                case 2; STmat = cat(1,statmat{2,ichh})'; % stimon Failed vs Succes
                    cols = colors(3:4,:);
                    tit = p.groups2comp{1,2}; leg = {p.groups2comp{2,1}, p.groups2comp{2,2}};
                case 3; s1 =statmat{1,ichh}; s2 = statmat{2,ichh};
                    STmat = cat(1,s1(1), s2(1)); % stimoff Failed vs stimon Failed
                    cols = colors([1 3],:);
                    tit = p.groups2comp{2,1}; leg = {p.groups2comp{1,1}, p.groups2comp{1,2}};
                case 4; s1 =statmat{1,ichh}; s2 = statmat{2,ichh};
                    STmat = cat(1,s1(2), s2(2));% stimoff Succes vs stimon Succes
                    cols = colors([2 4],:);
                    tit = p.groups2comp{2,2}; leg = {p.groups2comp{1,1}, p.groups2comp{1,2}};
            end;
        end
        
        pow = cellfun(@(x) squeeze(mean(x,1)),STmat,'UniformOutput',0);
        
        pmat = pcond{j,ichh};
        psig = any(pmat{1},1);
        
        subplot(compnr,compnr,j)
        for k = 1:length(pow)
            
            if strcmp(p.errorshade,'SE'); % standard error
                shad = (std(pow{k},[],2)) / sqrt(length(pow{k}));
            else strcmp(p.errorshade,'SD'); % standard deviation
                shad = std(pow{k},[],2);
            end
            he{k} = errorshade(new_times,mean(pow{k},2),shad,...
                'LineColor',cols(k,:),'ShadeColor',cols(k,:));hold on;
            
        end
        ylim(p.avg_clim);
        yL = p.avg_clim;
        scatter(new_times(psig), ones( 1,sum(psig) )*yL(1) ,'filled','k');
        %                 legend(arrayfun(@(x) [p.groups2comp{x} ', n=' num2str(patnrs(x))]  ,1:2,'UniformOutput',0))
        xlabel(['Time relative to ' p.event2comp ' (s)']);
        ylabel('Norm. pow');
        frnm = p.fr_name; frnm(strfind(frnm,'_'))= ' ';
        title({tit,frnm})
        legend([he{1}(2) he{2}(2)], leg);
        j = j+1;
    end
    try
        suptitle(chanlabels{ichh});
    catch
        disp('')
    end
    
    figdir = fullfile(resdir_statfigs,'pow_timeseries'); if ~isfolder(figdir); mkdir(figdir); end;
    fnm = fullfile(figdir,[ p.sessgroup '_' chanlabels{ichh} '_' num2str(p.alpha) '_' p.baseline_type]);
    
    saveas(fig,[fnm '.jpg'])
    saveas(fig,[fnm '.fig'])
    close(fig);
end
end



%--------------------------------------------------------------------------
function bargraphs_statmat(plotwin,plotfr,groups2comp_labels,groups2comp,event2comp,varargin)

prs = inputParser;

addRequired(prs,'plotwin',@isvector);
addRequired(prs,'plotfr',@isvector);
addRequired(prs,'groups2comp_labels',@iscell);
addRequired(prs,'groups2comp',@iscell);
addRequired(prs,'event2comp',@ischar);
addParameter(prs,'sessgroup','left',@ischar);

addParameter(prs,'topo_wins',[-0.5 0; 0 0.5],@ismatrix);
addParameter(prs,'stat','eeglab',@ischar);
addParameter(prs,'alpha',0.01,@isnumeric);
addParameter(prs,'mcorrect','fdr',@ischar);

addParameter(prs,'fr_name','delta',@ischar);
addParameter(prs,'freq_bands',[1 4],@isvector);

addParameter(prs,'rectime','postop',@ischar);
addParameter(prs,'rectype','EEG',@ischar);
addParameter(prs,'side','left',@ischar);
addParameter(prs,'csd',true,@islogical);
parse(prs,plotwin,plotfr,groups2comp_labels,groups2comp,event2comp,varargin{:})
p = prs.Results;



global rootdir figdir_pd

wnnm = num2str(p.topo_wins);
wnnm = wnnm(:)';

new_sr = 50;
stat_win = [p.topo_wins(1,1) p.topo_wins(end,end)];
stat_times = stat_win(1):1/new_sr:stat_win(2);

plot_times = plotwin(1):1/new_sr:plotwin(2);
plotinx = dsearchn(stat_times',plot_times');


load(fullfile(rootdir,'freq_components.mat'));
f_ind_st= intersect(find(f>=p.freq_bands(1)),find(f<=p.freq_bands(2)));
f_stat = f(f_ind_st);
f_ind_plot= intersect(find(f_stat>=p.plotfr(1)),find(f_stat<=p.plotfr(2)));


compnr = length(p.groups2comp_labels);
if compnr==1
    
    gr2c = [p.groups2comp{1} '_vs_' p.groups2comp{2}];
else
    
    gr2c = [p.groups2comp{1,1} '_vs_' p.groups2comp{1,2} '_' p.groups2comp{2,1} '_' p.groups2comp{2,2}];
end


resdir_statfigs = fullfile(figdir_pd,[p.rectime '_' p.rectype],'stats',p.sessgroup,...
    p.fr_name,p.event2comp,gr2c,['CSD' char(string(p.csd)) '_' p.stat '_Win' wnnm]);


if strcmp(p.rectype,'LFP')
    load(fullfile(resdir_statfigs,['statmatx_' p.side '_' num2str(p.alpha) '_' p.subregion 'STN_' p.baseline_type '.mat']));
else
    load(fullfile(resdir_statfigs,['statmatx_' p.side '_' num2str(p.alpha) '_' p.baseline_type '.mat']));
end

channr = length(chanlabels);
for ichh = 1:channr
    fig = figure;
    set(fig,'Visible','off')
    
    colors = [1 0 0; 0 1 0; 0 0 0; 0 0 1];
    if compnr==1
        stval = 1;
    else
        stval = 1:4;
    end
    
    j = 1;
    for st = stval
        
        if compnr==1
            STmat = statmat(:,ichh);
            leg = {p.groups2comp{1,1}, p.groups2comp{1,2}};
            tit = '';
            cols = colors(1:2,:);
        else
            switch st; case 1; STmat = cat(1,statmat{1,ichh})'; % stimoff Failed vs Succes
                tit = p.groups2comp{1,1}; leg = {p.groups2comp{2,1}, p.groups2comp{2,2}};
                case 2; STmat = cat(1,statmat{2,ichh})'; % stimon Failed vs Succes
                    tit = p.groups2comp{1,2}; leg = {p.groups2comp{2,1}, p.groups2comp{2,2}};
                case 3; s1 =statmat{1,ichh}; s2 = statmat{2,ichh};
                    STmat = cat(1,s1(1), s2(1)); % stimoff Failed vs stimon Failed
                    tit = p.groups2comp{2,1}; leg = {p.groups2comp{1,1}, p.groups2comp{1,2}};
                case 4; s1 =statmat{1,ichh}; s2 = statmat{2,ichh};
                    STmat = cat(1,s1(2), s2(2));% stimoff Succes vs stimon Succes
                    tit = p.groups2comp{2,2}; leg = {p.groups2comp{1,1}, p.groups2comp{1,2}};
            end;
        end
        
        pow = cellfun(@(x) squeeze(mean(x(f_ind_plot,plotinx,:),[1 2])),STmat,'UniformOutput',0);
        bp = cat(2,pow{:});
        
        
        subplot(compnr,compnr,j)
        
        boxplot(bp,leg); hold on;
        for kk = 1:size(bp,2);
            for kkk = 1:size(bp,1)
                scatter(kk+randi([-100 100],1)*0.0015,bp(kkk,kk),[],'k','filled'); hold on;
            end
        end
        set_my_boxplot(gca)
        yL = ylim;
        [pval, ~, ~] = ranksum(bp(:,1),bp(:,2));
        
        if pval<0.05; col = 'r'; else; col = 'k'; end;
        text(1,yL(2)*.9,['ranksum p=' num2str(pval)],'Color',col);
        
        
        [~, pval2, ~] = ttest2(bp(:,1),bp(:,2));
        if pval2<0.05; col = 'r'; else; col = 'k'; end;
        text(1,yL(2)*.75,['ttest p=' num2str(pval2)],'Color',col);
        
        title({[num2str(plotwin) ' s relative to ' p.event2comp ],[num2str(p.plotfr) 'Hz' ], tit});
        ylabel('Norm. pow');
        j = j+1;
    end
    
    suptitle(chanlabels{ichh});
    set(gcf,'Position',get(0,'Screensize'))
    
    figdir = fullfile(resdir_statfigs,'bar_stats',['WIN' num2str(p.plotwin) '_FR' num2str(plotfr)]); if ~isfolder(figdir); mkdir(figdir); end;
    fnm = fullfile(figdir,['WIN' num2str(p.plotwin) '_FR' num2str(plotfr) '_'  p.sessgroup '_' chanlabels{ichh} ]);
    fnm = [fnm '_' p.baseline_type];
    saveas(fig,[fnm '.jpg'])
    saveas(fig,[fnm '.fig'])
    close(fig);
end
end

%-----------------------------------------------------------------------------
function within_pat_stat_loop(sess2analyse,EventTypes,SubEventTypes,freq_bands,statwin,subregion,alpha,mcorrect,avg_clim,diff_clim,baseline_win,csd)

global rootdir figdir_pd
rectype = sess2analyse(1).rectype;
rectime = sess2analyse(1).rectime;



if strcmp(rectype,'EEG') && strcmp(rectime,'postop')
    load(fullfile(rootdir,'postop_EEG_chanlocs.mat'));
    chanlabels = {chanlocs.labels};
elseif strcmp(rectype,'EEG') &&  strcmp(rectime,'intraop')
    chanlabels = {'F4'};
elseif strcmp(rectype,'LFP')
    chanlabels = {'Ch1','Ch2','Ch3','Ch4','Ch5'};
end


for si = 1:length(sess2analyse)
    
    curr_resdir = sess2analyse(si).folder;
    patnm = sess2analyse(si).patient;
    tag =sess2analyse(si).tag;
    side =sess2analyse(si).side;
    
    fprintf('%s, %s...',patnm,side);
    
    for ei = 1:length(EventTypes)
        event = EventTypes{ei};
        subevs = SubEventTypes(ei,:);
        fprintf('%s...',event);
        
        resdir = fullfile(figdir_pd,[rectime '_' rectype],'stats','Within_pat',[event 'Failed_vs_Succesf'],...
            ['eeglab_Win' num2str(statwin)]);
        
        
        for ci = 1:length(chanlabels)
            act_chan = chanlabels{ci};
            
            fprintf('%s...',act_chan);
            
            within_pat_stat(event,subevs,patnm,side,tag,curr_resdir,act_chan,freq_bands,statwin,...
                rectype,subregion,resdir,alpha,mcorrect,avg_clim,diff_clim,baseline_win,csd)
        end
    end
    fprintf('\n');
    
end
end



%--------------------------------------------------------------------------
function within_pat_stat(event,subevs,patnm,side,tag,curr_resdir,act_chan,freq_bands,statwin,rectype,subregion,resdir,alpha,mcorrect,avg_clim,diff_clim,baseline_win,csd)

global rootdir


sr = 50;
eptime = -2:1/sr:2; eptime = eptime(1:end-1);

pltime = statwin(1):1/sr:statwin(2);
timinx = dsearchn(eptime',pltime');


load(fullfile(curr_resdir,'Evinxx.mat'));
load(fullfile(rootdir,'freq_components.mat'));
% Load TF data for one channel (freq x time) - contains whole recording

for sei = 1:2
    evty = subevs{sei};
    epodir = fullfile(curr_resdir, ['EventAVGs_CSD' char(string(csd))]);
    epoch_nm = [ event '_' evty '_EPOCHs.mat'];
    
    load(fullfile(epodir,epoch_nm))
    
    if ~strcmp(act_chan,'chanmean')
        evpow_TRN = EventEPOCHs.(act_chan).TrialNorm;
    else
        chanlabs = {'Ch1','Ch2','Ch3','Ch4','Ch5'};
        evpow_TRN_0 = cell(1,length(chanlabs));
        for ci = 1:length(chanlabs)
            if isfield(EventEPOCHs,chanlabs{ci})
                evpow_TRN_0{ci} = EventEPOCHs.(chanlabs{ci}).TrialNorm;
            end
        end
        evpow_TRN = nanmean(cat(4,evpow_TRN_0{:}),4);
        
    end
    evpow_subs{sei,1} = evpow_TRN(:,timinx,:)
    trnr{sei} = size(evpow_TRN,3);
end

if ~isempty(baseline_win)
    baslim_inx = dsearchn(pltime',baseline_win'); basinx = baslim_inx(1):baslim_inx(2);
    basdat = cat(3,evpow_subs{:});
    basavg = nanmean(basdat(:,basinx,:,:),[2 3]);
    bassd = std(basdat(:,basinx,:,:),[],[2 3],'omitnan');
    repbas =  repmat( basavg,[1 length(timinx) 1]);
    repbas_sd =  repmat( bassd,[1 length(timinx) 1]);
    notempty = ~cellfun(@isempty,evpow_subs);
    evpow_subs(notempty) = cellfun(@(x) (x-repbas./repbas_sd) ,evpow_subs(notempty),'UniformOutput',0);
end

for fri = 1:size(freq_bands,1)
    freqs = freq_bands(fri,:);
    
    
    resdirfr = [resdir '_FR' num2str(freqs)];
    if ~isfolder(resdirfr); mkdir(resdirfr); end;
    
    f_ind = find(f>freqs(1)&f<freqs(2));
    ff = f(f_ind);
    
    Lens = cellfun(@(x) size(x,3),evpow_subs);
    if ~any(Lens<2)
        evpow_subs2 = cellfun(@(x) x(f_ind,:,:),evpow_subs,'UniformOutput',0);
        
        %         [pcond, ~, ~, statscond, ~, ~] = std_stat(evpow_subs2,'condstats','on','groupstats','off','mode','eeglab',...
        %             'method','perm','mcorrect',mcorrect,'alpha',alpha);
        
        [F, df, pval]  = statcondfieldtrip(evpow_subs2, 'method','montecarlo','naccu', 1000,...
            'mcorrect',mcorrect,'alpha',alpha);
        pcond{1} = pval<alpha;
        statscond{1} = F;
        
        
    else
        fprintf('Not enough %s, %s %s',subevs{Lens<2},patnm,side);
        return;
    end
    
    
    
    fig = figure;
    set(fig,'Visible','off')
    for si = 1:2
        subplot(3,1,si);
        tf = mean(evpow_subs2{si,1},3);
        spectr_fig(tf,ff,sr,freqs(2),freqs(1),statwin, [], 0,avg_clim);
        title([subevs{si} '(n=' num2str(trnr{si}) ' trials)' ]);
        xL = xlim; xticks(round(linspace(xL(1),xL(end)-0.5,5))); xticklabels(arrayfun(@num2str,pltime(xticks),'UniformOutput',0));
        xlabel('Time (s)');
        
    end
    
    subplot(3,1,3)
    tf3 = statscond{1};
    spectr_fig(tf3,ff,sr,freqs(2),freqs(1),statwin, [], 0,diff_clim);
    hold on;
    contour(pcond{1},'Color','white'); title('Difference (F-S)')
    
    xL = xlim; xticks(round(linspace(xL(1),xL(end)-0.5,7))); xticklabels(arrayfun(@num2str,pltime(xticks),'UniformOutput',0));
    xlabel('Time (s)');
    
    suptitle({[patnm ' ' side],['Aligned to ' event],act_chan,[num2str(freqs) ' Hz']})
    set(fig,'Position',get(0,'Screensize'))
    
    fnm = [patnm '_' side '_' tag '_' act_chan  '_' num2str(alpha) '_MCORR' mcorrect 'BAS' num2str(baseline_win)];
    
    if strcmp(rectype,'LFP')
        saveas(fig,fullfile(resdirfr,[fnm '_' subregion 'STN.jpg']));
        saveas(fig,fullfile(resdirfr,[fnm '_' subregion 'STN.fig']));
%         saveas(fig,fullfile(resdirfr,[fnm '_' subregion 'STN.emf']));
        close(fig)
    else
        saveas(fig,fullfile(resdirfr,[fnm '.jpg']));
        saveas(fig,fullfile(resdirfr,[fnm '.fig']));
        %         saveas(fig,fullfile(resdirfr,[fnm '.emf']));
        close(fig)
    end
    
    %
    %     fig = figure;
    %
    %     tf3 = statscond{1};
    %     spectr_fig(tf3,ff,sr,freqs(2),freqs(1),statwin, [], 0,diff_clim);
    %     colormap(cool)
    %     hold on;
    %     contour(pcond{1},'Color','white'); title('Difference (F-S)')
    %
    %     xL = xlim; xticks(round(linspace(xL(1),xL(end)-0.5,7))); xticklabels(arrayfun(@num2str,pltime(xticks),'UniformOutput',0));
    %     xlabel('Time (s)');
    %
    %     suptitle({[patnm ' ' side],['Aligned to ' event],act_chan,[num2str(freqs) ' Hz']})
    %     set(fig,'Position',get(0,'Screensize'))
    %
    %     fnm = [patnm '_' side '_' tag '_' act_chan  '_' num2str(alpha) '_MCORR' mcorrect '_DIFFMAP'];
    %     if strcmp(rectype,'LFP')
    %         saveas(fig,fullfile(resdirfr,[fnm '_' subregion 'STN.jpg']));
    %         saveas(fig,fullfile(resdirfr,[fnm '_' subregion 'STN.fig'])); close(fig)
    %     else
    %         saveas(fig,fullfile(resdirfr,[fnm '.jpg']));
    %         saveas(fig,fullfile(resdirfr,[fnm '.fig'])); close(fig)
    %     end
    
    
end

end
