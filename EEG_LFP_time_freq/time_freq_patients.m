function time_freq_patients(sess2analyse,EventTypes,SubEventTypes,epoch_win,baseline_win,varargin);
%TIME_FREQ_PATIENTS     Time-frequency decomposition with wavelet 
%   TIME_FREQ_PATIENTS(sess2analyse,EventTypes,SubEventTypes,epoch_win,baseline_win,...)
%       -Applies wavelet transformation (see eegwavelet2.m) to preprocessed 
%       and epoched data (see preprocess_PD).
%       Data epochs are concatenated and zscored prior to transformation
%       (wavelet analysis is based on the code of Christopher Torrence and Gilbert P. Compo,1995-1998).
%       Complex wavelet coefficients and frequency vectors are saved in results 
%       folder of each patient in 'TFpows_blocks' subdirectory for each channel separately.
%       Wavelet coeficcients of all epochs are saved in 4 blocks (for faster loading of data).
%       
%       - Performs Full-epoch normalization and Baseline normalization for
%       sinlge epochs (separately). Performs permutation test with cluster-based correction across
%       epochs, within each patient, for each channel separately (the two type of normalized data is
%       tested separately). The resulting data are stored in structs.
%       Full-epoch normalized data is stored under the fieldname 'TrialNorm' field, 
%       baseline normalized data is stored under the fieldname 'BAS...'.
%       Normalized epochs are saved as '*EPOCHs.mat', -epoch averages
%       as '*AVGs.mat', -the result of the statistical tests as '*STAT.mat' files. 
%       Data resulting from LFP channel averages are saved in separate files (same file names, only extended with '_chan')
%       The file names contain the corresponding event and subevent label (ex: 'StimulusOn_FailedStopTrial_AVGs.mat').
%       The files are saved in 'EventAVGs_...' subfolder in each patient's result directory (SESS2ANALYSE.FOLDER). 
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
% Optional inputs (name-value pairs)
%   'stat_time'     1x2 vector, time window relative to event timestamp in sec, to perform statistics
% 
%   'csd'           true | false, if true, CSD transformed EEG data is used
%                   (relevant only for postoperative EEG data)
% 
%   'bipol'         true | false, if true, F4-F3 bipolar derivation of EEG data is used
%                   (relevant only for postoperative EEG data)
% 
%   'alphas'        numeric or vector, list of alpha threshold to consider when performing statistics
%                   ERSP masks (1 where stat. significant, 0 otherwise) are created using the the specified alpha threshold
% 
% See also: preprocess_PD, EEGWAVELET2

% Johanna Petra Szabó, Hangya Balázs, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu



prs = inputParser;
addRequired(prs,'sess2analyse',@isstruct);
addRequired(prs,'EventTypes',@iscell);
addRequired(prs,'SubEventTypes',@iscell);
addRequired(prs,'epoch_win',@isvector);
addRequired(prs,'baseline_win',@isvector);
addParameter(prs,'stat_time',[-1 1],@isvector);
addParameter(prs,'csd',true,@islogical);
addParameter(prs,'bipol',false,@islogical);
addParameter(prs,'choi',{},@iscell);
addParameter(prs,'alphas',NaN,@(x) isnumeric(x)||isvector(x));
parse(prs,sess2analyse,EventTypes,SubEventTypes,epoch_win,baseline_win,varargin{:});
p = prs.Results;



rectype = sess2analyse(1).rectype;


% Wavelet transformation
PD_wav(sess2analyse,p.choi,'csd',p.csd,'bipol',p.bipol); % wavelet of epochs(epochs concatenated)



% Full-epoch length single trial normalization + power averages across epochs (patient-by-patient)
save_patavgs_trialnorm(sess2analyse,EventTypes,SubEventTypes,epoch_win,'csd',...
    p.csd,'bipol',p.bipol,'alpha',p.alphas);

% Baseline single trial normalization + power averages across epochs (patient-by-patient)
save_patavgs_bascorr(sess2analyse,EventTypes,SubEventTypes,epoch_win,...
    baseline_win,'stat_win',p.stat_time,'csd',p.csd,'bipol',p.bipol,'measure','Pow','alpha',p.alphas);


if strcmp(rectype, 'LFP')
    save_patavgs_bascorrCHANMEAN(sess2analyse,EventTypes,SubEventTypes,epoch_win,...
        baseline_win,'stat_win',p.stat_time,'csd',p.csd,'alpha',p.alphas);
end


end



%--------------------------------------------------------------------------
function PD_wav(sess2analyse,choi,varargin)
%PD_WAV     Wavelet transform of all data
% Applies wavelet transformation (see eegwavelet2) to preprocessed epoched data (see preprocess_PD).
% Trial data are concatenated and zscored prior to transformation
% (wavelet analysis is based on the code of Christopher Torrence and Gilbert P. Compo,1995-1998).
% 
% Complex wavelet coefficients and frequency vectors are saved in results folder of each patient in 'TFpows_blocks' subdirectory.
% Epoched data is saved in 4 blocks (for faster loading of data).

% 
% Required inputs:
%   SESS2ANALYSE - struct containing all necessary information (name of patient, side
%     of experiment, tag of condition, session folder path) of sessions that
%     need to be analysed (see getdata2analyse)
%
%   CHOI - cell array of channels of interest; if empty, loops over all
%     channels
% 
% Optional input:
%   'csd'           true | false, if true, CSD transformed EEG data is used
%                   (relevant only for postoperative EEG data)
% 
%   'bipol'         true | false, if true, F4-F3 bipolar derivation of EEG data is used
%                   (relevant only for postoperative EEG data)

%
% See also: EEGWAVELET2


global rootdir


prs = inputParser;
addRequired(prs,'sess2analyse',@isstruct)
addRequired(prs,'choi',@iscell)
addParameter(prs,'csd',false,@islogical);
addParameter(prs,'bipol',false,@islogical);


parse(prs,sess2analyse,choi,varargin{:});
p = prs.Results;

%%
for snr = 1:length(sess2analyse)
    curr_resdir = sess2analyse(snr).folder;
    
    side = sess2analyse(snr).side;
    tag = sess2analyse(snr).tag;
    patnm = sess2analyse(snr).patient;
    
    fprintf('%s %s %s...\n',patnm, side, tag);
    
    try
        if p.csd
            EEG_ep = pop_loadset([curr_resdir filesep 'EEG_2plot_CSD.set']);
            tf_dir = fullfile(curr_resdir, 'TFpows_blocks_CSD');
            fastif(~isdir(tf_dir),mkdir(tf_dir),0);
        elseif p.bipol
            EEG_ep = pop_loadset([curr_resdir filesep 'EEG_2plot_bipol.set']);
            tf_dir = fullfile(curr_resdir, 'TFpows_blocks_bipol');
            fastif(~isdir(tf_dir),mkdir(tf_dir),0);
        else
            EEG_ep = pop_loadset([curr_resdir filesep 'EEG_2plot.set']);
            tf_dir = fullfile(curr_resdir, 'TFpows_blocks');
            fastif(~isdir(tf_dir),mkdir(tf_dir),0);
            
        end
    catch
        fprintf('EEG_ep fail: %s,%s,%s.\n',patnm,side,tag);
        continue
    end
    
    
    % Wavelet transform
    try
        [EEG_ep_choi] = concat_wav_pd(choi,EEG_ep,4,tf_dir);
        
    catch
        fprintf('Wavtransf fail: %s,%s,%s.\n',patnm,side,tag);
        continue
    end
end
end




%--------------------------------------------------------------------------
function save_patavgs_trialnorm(sess2analyse,EventTypes,SubEventTypes,epoch_win,varargin)

%%
dbstop if error

prs = inputParser;
addRequired(prs,'sess2analyse',@isstruct);
addRequired(prs,'EventTypes',@iscell);
addRequired(prs,'SubEventTypes',@iscell);
addRequired(prs,'epoch_win',@isvector);
addParameter(prs,'csd',false,@islogical);
addParameter(prs,'bipol',false,@islogical);
addParameter(prs,'measure','Pow',@ischar);
addParameter(prs,'chan','F4',@ischar);
addParameter(prs,'alpha',NaN,@isnumeric);
addParameter(prs,'mcorrect','fdr',@ischar); % fdr  | none
parse(prs,sess2analyse,EventTypes,SubEventTypes,epoch_win,varargin{:})
p = prs.Results;



global rootdir 


rectime = p.sess2analyse(1).rectime;
rectype = p.sess2analyse(1).rectype;

if ~strcmp(rectime,'postop')
    p.csd = false; p.bipol = false;
end

if strcmp(rectype,'EEG')
    chanlocs = [];
    load(fullfile(rootdir,'postop_EEG_chanlocs.mat'));
    if p.bipol
        chanlabels = {'F4-F3'};
    else
        chanlabels = {chanlocs.labels};
    end
    if strcmp(rectime,'intraop')
        finx = strcmp(p.chan,chanlabels);
        chanlocs = chanlocs(finx);
        chanlabels = {p.chan};
    end
else
    chanlabels = arrayfun(@(x) ['Ch' num2str(x)] ,1:5,'UniformOutput',0);
end



if p.csd
    tfdirnm = 'TFpows_blocks_CSD';
elseif p.bipol
    tfdirnm = 'TFpows_blocks_bipol';
else
    tfdirnm = 'TFpows_blocks';
end


Evinxx = [];

if isempty(p.SubEventTypes); subnr = 3; else; subnr = 1:3; end;

new_sr = 50;
new_time = epoch_win(1):1/new_sr:epoch_win(2);


%% Patient loop (in one group)
sessnr = length(sess2analyse);

for si = 1:sessnr
    
    curr_resdir = sess2analyse(si).folder;
    if p.bipol
        
        resdir = fullfile(curr_resdir, 'EventAVGs_bipol');
    else
        resdir = fullfile(curr_resdir, ['EventAVGs_CSD' char(string(p.csd))]);
    end
    if ~isfolder(resdir); mkdir(resdir); end;
    fprintf('%s...',sess2analyse(si).patient)
    
    
    % Channel loop
    chnr = length(chanlabels);
    for ich = 1:chnr
        act_chan = chanlabels{ich};
        fprintf('\n%s...',act_chan)
        
        
        
        
        % Load TF data for one channel (freq x time) - contains whole recording
        
        [epoch_pow, f] = load_TFblocks(curr_resdir,act_chan,p.measure,fullfile(sess2analyse(si).folder,tfdirnm),4);
        
        if isempty(epoch_pow)
            fprintf('Empty epoch_pow %s\n',curr_resdir)
            continue
        end
        
        
        
        %% EventTyepe loop
        for ei = 1:length(EventTypes)
            event = EventTypes{ei};
            fprintf('%s...',event)
            
            for sei = subnr
                if sei<3
                    
                    evty = SubEventTypes{ei,sei};
                    
                    fprintf('%s...',evty)
                else
                    evty = event;
                end
                
                % Find indeces of respective event type from whole recording
                evpow = [];
                
                
                
                
                evavg_nm = [ event '_' evty '_AVGs.mat'];
                epoch_nm = [ event '_' evty '_EPOCHs.mat'];
                evSTAT_nm = [ event '_' evty '_STATs.mat'];
               
                
                
                if exist(fullfile(resdir,evavg_nm))==2
                    load(fullfile(resdir,evavg_nm))
                else
                    EventAVGs = struct;
                end
                
                if exist(fullfile(resdir,evSTAT_nm))==2
                    load(fullfile(resdir,evSTAT_nm))
                else
                    EventSTAT = struct;
                end
                
                  
                if exist(fullfile(resdir,epoch_nm))==2
                    load(fullfile(resdir,epoch_nm))
                else
                    EventEPOCHs = struct;
                end
                
                
                
                
                
                load(fullfile(curr_resdir,'Evinxx.mat'));
                
                
                
                if ~strcmp(event,'StopSignal') && ismember(evty,{'FailedStopTrial','SuccesfulStopTrial'});
                    
                    [~,evinx] = StimOn_stoppart_evinx(Evinxx,event,evty);
                else
                    evinx = Evinxx.(event).(evty).epoch_index;
                end
                
                if isempty(evinx)
                    fprintf('No %s trial.\n',evty);
                    continue
                end
                
                try
                    evpow = epoch_pow(:,:,evinx); 
                catch
                    fprintf('evinx mismatch %s \n', curr_resdir)
                end
                
                
                
                
                
                % NORMALIZATION
              
                    
                    if p.bipol
                        act_chan(strfind(act_chan, '-')) = '_';
                    end
                    % Full-epoch length single-trial norm. (TrialNorm)
                    M = mean(evpow,2);
                    SD = std(evpow,[],2);
                    Mrep = repmat(M,[1 size(evpow,2)]);
                    SDrep = repmat(SD,[1 size(evpow,2)]);
                    
                    evpowTRN = (evpow - Mrep)./SDrep;
                    
                    TRN_AVG = nanmean(evpowTRN,3); % TrialNorm normalized trials averaged
                    EventEPOCHs.(act_chan).TrialNorm =  evpowTRN;
                    EventAVGs.(act_chan).TRN_AVG =  TRN_AVG;
                    
                    formula = 'mean(arg1,3);';
                    
                    
                
                
                % STATISTICS
                
                if ~isnan(p.alpha)
                    hold on;
                    [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(evpowTRN,f,p.alpha,1000,false,p.mcorrect,formula);
                    EventSTAT.(act_chan).TrialNorm.p_ersp =  exactp_ersp;
                    EventSTAT.(act_chan).TrialNorm.mask_ersp =  maskersp;
                    EventSTAT.(act_chan).TrialNorm.alphafdr =  alphafdr;
                    EventSTAT.(act_chan).TrialNorm.alpha =  p.alpha;
                    EventSTAT.(act_chan).TrialNorm.trialnr =  size(evpowTRN,3);
                    
                end
                
                
                
                
                
                save(fullfile(resdir,evavg_nm),'EventAVGs');
                save(fullfile(resdir,evSTAT_nm),'EventSTAT');
                save(fullfile(resdir,epoch_nm),'EventEPOCHs');
                
                
                
            end
        end
    end
    
    %     end
    
end
end





%--------------------------------------------------------------------------
function save_patavgs_bascorr(sess2analyse,EventTypes,SubEventTypes,epoch_win,baseline_win,varargin)

%%
dbstop if error

prs = inputParser;
addRequired(prs,'sess2analyse',@isstruct);
addRequired(prs,'EventTypes',@iscell);
addRequired(prs,'SubEventTypes',@iscell);
addRequired(prs,'epoch_win',@isvector);
addRequired(prs,'baseline_win',@(x) isvector(x)| isnumeric(x));
addParameter(prs,'stat_win',[],@(x) isvector(x)| isnumeric(x));
addParameter(prs,'csd',false,@islogical);
addParameter(prs,'bipol',false,@islogical);
addParameter(prs,'measure','Pow',@ischar);
addParameter(prs,'chan','F4',@ischar);
addParameter(prs,'alpha',NaN,@isnumeric);
addParameter(prs,'mcorrect','fdr',@ischar); % fdr  | none
parse(prs,sess2analyse,EventTypes,SubEventTypes,epoch_win,baseline_win,varargin{:})
p = prs.Results;



global rootdir filesdir


rectime = p.sess2analyse(1).rectime;
rectype = p.sess2analyse(1).rectype;
if ~strcmp(rectime,'postop')
    p.csd = false; p.bipol = false;
end


Evinxx = [];

if isempty(p.SubEventTypes); subnr = 3; else; subnr = 1:3; end;

new_sr = 50;
new_time = epoch_win(1):1/new_sr:epoch_win(2); new_time = new_time(1:end-1);
load(fullfile(rootdir,'freq_components.mat'));

%% Patient loop (in one group)
sessnr = length(sess2analyse);

for si = 1:sessnr
    
    
    curr_resdir = sess2analyse(si).folder;
    if p.bipol
        
        resdir = fullfile(curr_resdir, 'EventAVGs_bipol'); 
    else
        resdir = fullfile(curr_resdir, ['EventAVGs_CSD' char(string(p.csd))]); 
    end
    if ~isfolder(resdir); mkdir(resdir); end;
    fprintf('%s...',sess2analyse(si).patient )
    
    
    %% EventTyepe loop
    for ei = 1:length(EventTypes)
        event = EventTypes{ei};
        
        fprintf('%s...',event )
        
        for sei = subnr
            if sei<3
                
                evty = SubEventTypes{ei,sei};
                
                fprintf('%s...',evty )
            else
                evty = event;
            end
            
            
            
            
            
           
            evavg_nm = [ event '_' evty '_AVGs.mat'];
            epoch_nm = [ event '_' evty '_EPOCHs.mat'];
            evSTAT_nm = [ event '_' evty '_STATs.mat'];
            
            try
                load(fullfile(resdir,evavg_nm))
                load(fullfile(resdir,epoch_nm))
            catch
                fprintf('No EventAVGs.mat\n' );
                continue;
            end
            load(fullfile(resdir,evSTAT_nm))
            
            chanlabs = fieldnames(EventAVGs);
            
            for ci = 1:length(chanlabs)
                
                act_chan = chanlabs{ci};
                
                if p.bipol
                    act_chan(strfind(act_chan, '-')) = '_';
                end
                
                try
                    evpowTRN = EventEPOCHs.(act_chan).TrialNorm;
                catch
                    continue;
                end
                fprintf('\n%s...',chanlabs{ci} )
                
                
                
                % NORMALIZATION
               
                    
                    baslims = dsearchn(new_time',p.baseline_win'); basinx = baslims(1):baslims(2);
                    statlims = dsearchn(new_time',p.stat_win');
                    
                    %                     B = repmat(nanmean(evpowTRN(:,basinx,:),[2]),[1 200 1]);
                    %                     Bsd = repmat(std(evpowTRN(:,basinx,:),[],[2],'omitnan'),[1 200 1]);
                    %                     evpowB = (evpowTRN - B) ./ Bsd; % Baseline correction for each trial for statistics
                    %                     B_AVG2 = nanmean(evpowB,3);
                    
                    B2 = repmat(nanmean(evpowTRN(:,basinx,:),[2 3]),[1 200]);
                    Bsd2 = repmat(std(evpowTRN(:,basinx,:),[],[2 3],'omitnan'),[1 200]);
                    %
                    B_AVG2 = (nanmean(evpowTRN,3)- B2)./ Bsd2; % Baseline correction for trial average for visualization
                    %                     B_AVG2 = nanmean(evpowTRN,3)./ B2;
                    
                    basnm = ['BASLims' num2str(baslims(1)) '_' num2str(baslims(2))];
                    EventAVGs.(act_chan).(basnm) =  B_AVG2(:,baslims(1):end);
                    
                    
                    formula = 'mean(arg1,3);';
                    
                
                
                
                
                % STATISTICS
                
                if ~isnan(p.alpha) && ~isempty(p.baseline_win)
                    hold on;
                    [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(evpowTRN(:,baslims(1):statlims(2),:),f,p.alpha,1000,false,p.mcorrect,formula,1:length(basinx));
               
                    
                    EventSTAT.(act_chan).(basnm).p_ersp =  exactp_ersp;
                    EventSTAT.(act_chan).(basnm).mask_ersp =  maskersp;
                    EventSTAT.(act_chan).(basnm).alphafdr =  alphafdr;
                    EventSTAT.(act_chan).(basnm).alpha =  p.alpha;
                    EventSTAT.(act_chan).(basnm).trialnr =  size(evpowTRN,3);
                end
                
                
                
                
                
                save(fullfile(resdir,evavg_nm),'EventAVGs');
                save(fullfile(resdir,evSTAT_nm),'EventSTAT');
                
                
                
            end
        end
    end
    
    %     end
    
end
end




%--------------------------------------------------------------------------
function save_patavgs_bascorrCHANMEAN(sess2analyse,EventTypes,SubEventTypes,epoch_win,baseline_win,varargin)

%%
dbstop if error

prs = inputParser;
addRequired(prs,'sess2analyse',@isstruct);
addRequired(prs,'EventTypes',@iscell);
addRequired(prs,'SubEventTypes',@iscell);
addRequired(prs,'epoch_win',@isvector);
addRequired(prs,'baseline_win',@(x) isvector(x)| isnumeric(x));
addParameter(prs,'stat_win',[],@(x) isvector(x)| isnumeric(x));
addParameter(prs,'csd',false,@islogical);
addParameter(prs,'alpha',NaN,@isnumeric);
addParameter(prs,'mcorrect','fdr',@ischar); % fdr  | none
parse(prs,sess2analyse,EventTypes,SubEventTypes,epoch_win,baseline_win,varargin{:})
p = prs.Results;



global rootdir filesdir





if isempty(p.SubEventTypes); subnr = 3; else; subnr = 1:3; end;

new_sr = 50;
new_time = epoch_win(1):1/new_sr:epoch_win(2); new_time = new_time(1:end-1);
load(fullfile(rootdir,'freq_components.mat'));

chanlabs = {'Ch1','Ch2','Ch3','Ch4','Ch5'};
%% Patient loop (in one group)
sessnr = length(sess2analyse);

for si = 1:sessnr
    
    
    curr_resdir = sess2analyse(si).folder;
  
        resdir = fullfile(curr_resdir, ['EventAVGs_CSD' char(string(p.csd))]);
   
    if ~isfolder(resdir); mkdir(resdir); end;
    fprintf('%s...',sess2analyse(si).patient )
    
    
    %% EventTyepe loop
    for ei = 1:length(EventTypes)
        event = EventTypes{ei};
        
        fprintf('%s...',event )
        
        for sei = subnr
            if sei<3
                
                evty = SubEventTypes{ei,sei};
                
                fprintf('%s...',evty )
            else
                evty = event;
            end
            
            
            epoch_nm = [ event '_' evty '_EPOCHs.mat'];
            evavg_nm = [ event '_' evty '_chan_AVGs.mat'];
            evSTAT_nm = [ event '_' evty '_chan_STATs.mat'];
            
            try
                load(fullfile(resdir,epoch_nm))
            catch
                fprintf('No EventEPOCH.mat\n' );
                continue;
            end
            
            if exist(fullfile(resdir,evavg_nm))==2
                load(fullfile(resdir,evavg_nm))
            else
                EventAVGs = struct;
            end
            
            if exist(fullfile(resdir,evSTAT_nm))==2
                load(fullfile(resdir,evSTAT_nm))
            else
                EventSTAT = struct;
            end
            
            
            evpowTRN_0 = cell(1,length(chanlabs));
            for ci = 1:length(chanlabs)
                
                act_chan = chanlabs{ci};
                
                                
                try
                    evpowTRN_0{ci} = EventEPOCHs.(act_chan).TrialNorm;
                catch
                    fprintf('no')
                    continue;
                end
                fprintf('\n%s...',chanlabs{ci} )
            end
            
            if isempty(evpowTRN_0)
                fprintf('No data %s %s\n',event,evty)
                continue
            end
            
            % CHANNEL MEAN
            evpowTRN= nanmean( cat(4,evpowTRN_0{:}) ,4);
            
            TRN_AVG = nanmean(evpowTRN,3);
            EventAVGs.chanmean.TRN_AVG =  TRN_AVG;
             
            % NORMALIZATION
            baslims = dsearchn(new_time',p.baseline_win'); basinx = baslims(1):baslims(2);
            statlims = dsearchn(new_time',p.stat_win');
            
            %                     B = repmat(nanmean(evpowTRN(:,basinx,:),[2]),[1 200 1]);
            %                     Bsd = repmat(std(evpowTRN(:,basinx,:),[],[2],'omitnan'),[1 200 1]);
            %                     evpowB = (evpowTRN - B) ./ Bsd; % Baseline correction for each trial for statistics
            %                     B_AVG2 = nanmean(evpowB,3);
            
            B2 = repmat(nanmean(evpowTRN(:,basinx,:),[2 3]),[1 200]);
            Bsd2 = repmat(std(evpowTRN(:,basinx,:),[],[2 3],'omitnan'),[1 200]);
            %
            B_AVG2 = (nanmean(evpowTRN,3)- B2)./ Bsd2; % Baseline correction for trial average for visualization

            
            basnm = ['BASLims' num2str(baslims(1)) '_' num2str(baslims(2))];
            EventAVGs.chanmean.(basnm) =  B_AVG2(:,baslims(1):end);
            
            
            formula = 'mean(arg1,3);';
            
            
            
            % STATISTICS
                
                
            if ~isnan(p.alpha) && ~isempty(p.baseline_win)
                hold on;
                [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(evpowTRN(:,baslims(1):statlims(2),:),f,p.alpha,1000,false,p.mcorrect,formula,1:length(basinx));
                
                
                EventSTAT.chanmean.(basnm).p_ersp =  exactp_ersp;
                EventSTAT.chanmean.(basnm).mask_ersp =  maskersp;
                EventSTAT.chanmean.(basnm).alphafdr =  alphafdr;
                EventSTAT.chanmean.(basnm).alpha =  p.alpha;
                EventSTAT.chanmean.(basnm).trialnr =  size(evpowTRN,3);
            end
            
            
            
            
            
            save(fullfile(resdir,evavg_nm),'EventAVGs');
            save(fullfile(resdir,evSTAT_nm),'EventSTAT');
            
            
            
        end
    end
end

%     end

end


%--------------------------------------------------------------------------
function [EEG_ep_choi] = concat_wav_pd(choi,EEG_ep,blocknr,tf_dir)
%CONCAT_WAV_PD Wavelet transform of a single EEG data


dbstop if error
fprintf('Wavelet transformation...\n')


if isempty(choi)
    try
    choi = {EEG_ep.chanlocs(:).labels};
    catch
        for i = 1:EEG_ep.nbchan
            EEG_ep.chanlocs(i).labels = ['Ch' num2str(i)];
        end
        choi = {EEG_ep.chanlocs(:).labels};
    end
end

EEG_ep_choi = pop_select(EEG_ep,'channel',choi);

blocklen = ceil(size(EEG_ep_choi.data,3)/blocknr);

ds = floor(linspace(1,size(EEG_ep.data,2),200));

c_min = [];
c_max = [];
for oc = 1:EEG_ep_choi.nbchan
    fprintf('Channel');
    act_chan = EEG_ep_choi.chanlocs(oc).labels;
    fprintf('%s ',act_chan);
    
    if exist([tf_dir filesep 'epoch_wav_ch' act_chan '_4.mat'])~=2
        signconcat =  reshape(EEG_ep_choi.data(oc,:,:),1,[]);
        
        Zsign = zscore(signconcat);
        
        % eegwavelet function is modified to give complex wavelet coefficients as output instead of power and phase values (due to memory shortage)
        [waveconcat,f] = eegwavelet2(Zsign,EEG_ep_choi.srate); 
        f = f(f>0.5);
        
        st = 1;
        ep_length = size(EEG_ep_choi.data,2);
        epnr = size(EEG_ep_choi.data,3);
        epoch_wav = nan(length(f),ep_length,epnr);
        for eo = 1:epnr
            try
                epoch_wav(:,:,eo) = waveconcat(:,st:st+ep_length-1);
            catch
                fprintf('error');
            end
            st = st+ep_length;
        end
        
        epoch_wavDS = epoch_wav(:,ds,:);
        
        % Save in blocks
        st2 = 1;
        for bi = 1:blocknr
            en = st2+blocklen-1;
            epoch_wav_bl = epoch_wavDS(:,:,st2:min(en,size(epoch_wav,3)));
            save(fullfile(tf_dir,['epoch_wav_ch' act_chan '_' num2str(bi) '.mat']),'epoch_wav_bl');
            
            st2 = st2+blocklen; epoch_wav_bl = [];
        end
        
        save(fullfile(tf_dir,['epoch_f_ch' act_chan '.mat']),'f');
        
    else
        
        fprintf('TF blocks done\n');
    end
    
end
end
