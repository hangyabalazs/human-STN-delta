function preprocess_PD(sess2analyse,epoch_win,EventTypes,SubEventTypes,preproc,baseline_win)
%PREPROCESS_PD preprocesses EEG/LFP data for further analysis
% PREPROCESS_PD(...) preprocesses raw, continuos EEG/LFP data, included in SESS2ANALYSE,
% saves filtered and cleaned data in EEGLAB (.set) format.
% This function requires the following toolboxes:
%       EEGLAB (Delorme A & Makeig S, 2004, 10.1016/j.jneumeth.2003.10.009)
%           + FileIO, ICLabel, Fieltrip-lite, firfilt, CleanLine plugins
%       DBSFILT (by Guillaume Lio,2012, https://github.com/guillaumelio/DBSFILT)
%       CSD (by Jürgen Kayser, 2009, doi:10.1016/j.clinph.2005.08.034)
%       FastICA (https://research.ics.aalto.fi/ica/fastica/)
%
% Preprocessing steps:
% 1. Load raw data, import channel data, convert into EEGLAB dataset, save
%              raw data is .set format.
% 2. Clean DBS artifacts (only if 'stimon' condition):
%              DBSFILT GUI (by Guillaume Lio,2012, https://github.com/guillaumelio/DBSFILT)
%              Applies frequency-domain Hampel filtering to clean DBS induced
%              artifacts.
% 3. Import behavioral data into EEG dataset (according to TrialEvents.mat struct
%              synchronized with EEG previously)
% 4. Resample: 250 Hz
% 5. Filter:   Lowpass filter with 100 Hz cutoff freq.
%              CleanLine filter (by Tim Mullen, 2011, https://github.com/sccn/cleanline)
%              Removes 50 Hz linenoise (uses sliding window to adaptively estimate sine wave amplitude to subtract)
%              If not sufficiently effective (has to be approved manually),
%              45 - 55 Hz notch filter is applied (pop_eegfiltnew, EEGLAB plugin).
%              Two types of highpass filter (resulting in two datasets):
%                 - 0.5 Hz cutoff freq - for data used for further analysis
%                 - 2 Hz cutoff freq - for ICA
%              Filtered datasets are saved in patient's result directory
%              (curr_resdir in sess2analyse struct)
% 6. Data epoching (only if preproc = 'epochs' !): data is epoched into trials,
%              according to behavioral events (EVENTTYPES, SUBEVENTTYPES)
%              and time window defined by EPOCH_WIN. Trial indeces are also saved
%              separately in Evinxx.mat file.
% 7. Baseline subtraction: baseline (BASELINE_WIN) is subtracted from each epoch, each channel.
% 8. Artifact rejection: - bad data portions/trials are removed manually from continuous/epoched EEG (pop_eegplot),
%                             -- 'continu' preprocessing (PREPROC): boundaries (index & duration of rejected data) are inserted to EEG.event field
%                             new EEG.event structure is saved in events_with_boundaries.mat file
%                             -- 'epochs' preprocessing (PREPROC): indeces of rejected trials are saved into rejected_epochs.mat file
%                        - (bad channels are removed and interpolated)
%                        - Independent Component Analysis (FastICA, https://research.ics.aalto.fi/ica/fastica/) 
%                        is performed following bad data/trial rejection
%                        using 2Hz highpass filtered data
%                        - Resulting components have to be reviewd visually.
%                        Then components corresponding to artifacts related to
%                        blinks/facial movements/bad channels are selected.
%                        Selected components to reject are also saved separately
%                        to 'gcompreject_continu.mat'
%                        - Selected components are subtracted from the 0.5 Hz
%                        filtered data.
% 9. Save final data structure: EEG_continu.set/EEG_2plot.set is used for further analysis.
%                               Channel data also saved into a separate file:
%                               [rectype '_EEG_chanlocs.mat']
% 10. Current Source Density transormation is applied to postoperative EEG data:
%                        CSD toolbox by Jürgen Kayser, 2009 (doi:10.1016/j.clinph.2005.08.034)
% 11. Converts final postop EEG data to bipolar montage and saves data of F4-F3 derivation.
%
% Required inputs:
%     SESS2ANALYSE      struct containing all necessary information (name of patient, side
%                       of experiment, tag of condition, session folder path) of
%                       session data that need to be analysed (see getdata2analyse)
%
%     EPOCH_WIN         1x2 vector, time window relative to event timestamp in sec, for data epoching (ex: [-2 2])
%
%     EVENTTYPES        1xN cell array of event labels, ex: {'StimulusOn','StopSignal','KeyPress1','Feedback'};
%
%     SUBEVENTTYPES     Nx2 cell array of partition ("subevent") labels, each row
%                       corresponds to an event label, columns to partitions
%                       {'CueStim','StopStim';'FailedStopTrial','SuccesfulStopTrial';'CueResponse','StopResponse';'Correct','Error';};
%
%     PREPROC           'continu'| 'epochs' - preprocess and save eeg in coninuous/
%                       epoched form
%
%     BASELINE_WIN       1x2 vector, time window relative to event timestamp in sec, for baseline correction
%
% See also:

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

%}

global rootdir ALLEEG CURRENTSET EEG ALLCOM


rectime = sess2analyse(1).rectime;
for snr = 1:length(sess2analyse)
    %%
    % Preallocate eeglab variables
    EEG =[]; ALLEEG = []; ALLCOM = {};  CURRENTSET = 0;
    
    
    curr_resdir = sess2analyse(snr).folder;
    side = sess2analyse(snr).side;
    condition = sess2analyse(snr).tag;
    patnm = sess2analyse(snr).patient;
    currsess = sess2analyse(snr).sessfolder;
    
    rectype = sess2analyse(snr).rectype;
    fprintf('%s %s %s...\n',patnm, side, condition);
    
    
    
    % Load raw data + filter/ Load filtered data
    
    [EEG_filt1 EEG_filt2 iste] = load_filter_data(currsess,patnm,side,condition,...
        curr_resdir,EventTypes,SubEventTypes,rectype, rectime);
    
    
    if iste==0 % if is TrialEvents.mat (if synchronization was successful)
        continue
    end
    
    if strcmp(preproc,'epochs')
        % Prepare epochs
        [EEG_filt1,EEG_filt2] =  prep_epochs(EEG_filt1, EEG_filt2, epoch_win);
        eeg_setnm = 'EEG_2plot.set';
        
%         
%     elseif strcmp(preproc,'continu')
%         eeg_setnm = 'EEG_continu.set';
    end
    
    
    %     Reject bad data + apply ICA
    
    if exist(fullfile(curr_resdir,eeg_setnm))~=2
        [EEG] = reject_bad_ica(EEG_filt1,EEG_filt2, curr_resdir);
        
        
        if strcmp(preproc,'epochs')
            % Save event indeces
            save_evinxx(EEG,EventTypes,SubEventTypes,curr_resdir,true);
            
            
            % Subtract baseline
            EEG = substr_bas(EEG,baseline_win);
        end
        
        
        
        % Save EEG 2 plot
        setnm = [curr_resdir filesep eeg_setnm];
        pop_saveset(EEG,setnm);
        chanlocs = EEG.chanlocs;
        save(fullfile(rootdir,[rectype '_EEG_chanlocs.mat']),'chanlocs')
    end
    
    
end

%% CSD transformation
if strcmp(rectime,'postop')
    EEG_CSD_ft(sess2analyse)
end



%% Re-reference 2 bipolar montage (F4-F3) to match intraop recording
if strcmp(rectime,'postop')
    reref2bipol(sess2analyse)
end



end

%--------------------------------------------------------------------------
function EEG = load_raw_data(currsess,rectype,rectime,patnm,side,condition,curr_resdir)

global  ALLEEG CURRENTSET EEG

switch rectype
    case 'EEG'
        
        switch rectime
            case 'postop'
                EEG = load_postopeeg(currsess,side,condition,curr_resdir);
                
            case 'intraop'
                EEG =  load_intraopeeg(currsess,patnm,side);
        end
        
    case 'LFP'
        EEG =  load_intraoplfp(currsess,patnm,side);
end

EEG.eegtype = [rectype '_' rectime];
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
end


%-------------------------------------------------------------------------
function EEG = load_postopeeg(currsess,side,tag,curr_resdir)
dbstop if error


EEG_ep1 = []; EEG_ep2 = []; iste = 1;

all_eegfile = dir([currsess filesep '*.eeg']);
currfile = find_filetag(all_eegfile,currsess,tag);

% Import data
%-------------
dbfilt = 0;
if ~isempty(dir([curr_resdir filesep 'EEG_raw*' '.set']))
    
    try
        %             EEG = pop_loadset([curr_resdir filesep 'EEG_raw_TemporalFiltered_DBSfiltered.set']);
        EEG = pop_loadset([currsess filesep tag '_' side filesep 'EEG_raw_TemporalFiltered_DBSfiltered.set']);
        fprintf('Loading DBS filtered data...\n');
        dbfilt = 1;
        
        
    catch
        %             EEG = pop_loadset([curr_resdir filesep 'EEG_raw.set']);
        EEG = pop_loadset([currsess filesep tag '_' side filesep 'EEG_raw.set']);
        fprintf('Loading raw data...\n');;
    end
    
else
    % Load .eeg file
    EEG = pop_fileio([currsess filesep currfile]);
    
    % Channel locations
    %--------------------
    EEG = pop_chanedit(EEG);
    
    %         Bad channel rejection
    %                 fig = figure;
    %                 pop_spectopo(EEG,1,[1,size(EEG.data,2)],'EEG','freqrange',[2 50],'title','Reject bad channel');
    %
    %         rejch = input('Press 1 for bad channel rejection, press any key otherwise.\n');
    %
    %         if rejch==1
    %             ch2rem = input('Nr. of channels to remove (for ex: [1,5,13])\n');
    ch2rem = [17,22,41,46];
    EEG = pop_select(EEG,'nochannel',ch2rem);
    %             EEG = pop_interp(EEG,ch2rem, 'spherical'); % interpolate
    %         end
    %                 close(fig)
    
    % Save file in eeglab dataset format
    setnm = [curr_resdir filesep 'EEG_raw.set'];
    pop_saveset(EEG,setnm);
end

if strcmp(tag,'stimon') && dbfilt==0
    diary on
    db = DBSFILT;
    keyboard
    diary off
    close(db)
    EEG = []; EEG = pop_loadset([curr_resdir filesep 'EEG_raw_TemporalFiltered_DBSfiltered.set']);
end
end



%--------------------------------------------------------------------------
function  EEG =  load_intraopeeg(currsess,patnm,side);

global rootdir
load(fullfile(rootdir,'sessioninfos_EEG_intraop.mat'));
patrow = find(strcmp(sessioninfos,patnm));

try
    eegfile = dir(fullfile(currsess,['*EMG*' sessioninfos{patrow,3} '*.mat']));
    eegdat = load(fullfile(eegfile.folder,eegfile.name));
catch
    fprintf('No EMG mat file %s\n',currsess);
    EEG = [];
    return
end

close(gcf)
EEG = pop_importdata('dataformat','array','nbchan',1,'data',eegdat.Data',...
    'setname', [patnm '_' side '_frontaleeg' ] ,'srate',eegdat.SampFreq,...
    'subject',patnm,'pnts',length(eegdat.t),'xmin',0);
% EEG.chanlocs.labels = sessioninfos{patrow,4};
% EEG.ref = sessioninfos{patrow,5};
EEG.chanlocs.labels = sessioninfos{patrow,5};
EEG.ref = sessioninfos{patrow,4};
end


%--------------------------------------------------------------------------
function [EEG_filt1 EEG_filt2 iste] = load_filter_data(currsess,patnm,side,condition,curr_resdir,EventTypes,SubEventTypes,rectype,rectime)



EEG_ep1 = []; EEG_ep2 = []; iste = 1;
cd(curr_resdir)

if isempty(dir([curr_resdir filesep 'EEG_*_2HP.set'])) && isempty(dir([curr_resdir filesep 'EEG_*_05HP.set']))
    
    %Load raw data
    %---------------
    
    EEG = load_raw_data(currsess,rectype,rectime,patnm,side,condition,curr_resdir);
    global  ALLEEG CURRENTSET EEG
    
    
    % Label events
    %--------------
    try
        EEG = behav_events(EEG,EventTypes,SubEventTypes,currsess,condition);
        iste = 1;
    catch
        fprintf('No TE file\n');
        iste = 0; EEG_filt1 = []; EEG_filt2 = [];
        return
    end
    
    % Downsample
    %-------------
    new_sr = 250;
    EEG = pop_resample(EEG,new_sr);
    EEG.srate = new_sr;
    
    % Filter
    %------------
    
    % Lowpass filter for stimoff
    if strcmp(condition,'stimoff')
        hicutoff = 100; % higher edge of passband: 100 Hz
        [EEG, com, b] = pop_eegfiltnew(EEG,[],hicutoff);
    end
    
    % Linenoise removal
    %%
    [EEG Lnfilt] = rem_line_noise(EEG);
    %%
    
    % High-pass filter: dataset 1 for analyses, dataset 2 for ICA
    
    locutoff = 0.5; % lower edge of passband: 0.5 Hz
    [EEG_filt1, com, b] = pop_eegfiltnew(EEG,locutoff,[]);
    
    
    locutoff = 2; % lower edge of passband: 2 Hz  (ICA might be biased by low freqs)
    [EEG_filt2, com, b] = pop_eegfiltnew(EEG,locutoff,[]);
    
    
    
    % Reject data chunk if there is long-lasting noise
    
    %     [EEG_filt1,EEG_filt2] = rej_badtrials(EEG_filt1,EEG_filt2,curr_resdir);
    
    
    % Save
    setnm = [curr_resdir filesep 'EEG_filt' Lnfilt '_05HP.set'];
    pop_saveset(EEG_filt1,setnm);
    
    setnm = [curr_resdir filesep 'EEG_filt' Lnfilt '_2HP.set'];
    pop_saveset(EEG_filt2,setnm);
    
    
else
    
    eegf2= dir([curr_resdir filesep 'EEG_*_2HP.set']);
    eegf05= dir([curr_resdir filesep 'EEG_*_05HP.set']);
    EEG_filt1 = pop_loadset(eegf05(1).name);
    EEG_filt2 = pop_loadset(eegf2(1).name);
    iste = 1;
end
end

%--------------------------------------------------------------------------
function reref2bipol(sess2analyse)
%REREF2BIPOL Re-reference data to bipolar montage (only F4-F3 derivation)
% REREF2BIPOL(sess2analyse)  Saves re-referenced EEG data as EEG_2plot_bipol.set

% Required input: sess2analyse: structure containing data to analyse (see getdata2analyse)
for snr = 1:length(sess2analyse)
    curr_resdir = sess2analyse(snr).folder;
    
    try
        EEG = pop_loadset(fullfile(curr_resdir,'EEG_2plot.set'));
    catch
        fprintf('NO EEG %s\n',curr_resdir);
        continue;
    end
    EEG = pop_select(EEG,'channel',{'F4','F3'});
    EEGbip = pop_reref(EEG,{'F3'});
    EEGbip.chanlocs(1).labels = 'F4-F3';
    
    pop_saveset(EEGbip,fullfile(curr_resdir,'EEG_2plot_bipol.set'))
    
    EEG = [];
    
end
end


%--------------------------------------------------------------------------
function EEG_CSD_ft(s2a)
%EEG_CSD_FT     Applies CSD tranformation on EEG data
%   EEG_CSD_ft(s2a) Applies CSD tranformation on EEG data defined
%   in S2A (see getdata2analyse), using CSD toolbox by Jürgen Kayser.
%       Kayser, J., Tenke, C.E. (2006a). doi:10.1016/j.clinph.2005.08.034
%       Kayser, J. (2009).Current source density (CSD) interpolation using spherical splines -
%               CSD Toolbox (Version 1.1) [http://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox].
%               New York State Psychiatric Institute: Division of Cognitive Neuroscience.
%
%   Transformed data is saved as EEG_2plot_CSD.set in the patient's result
%   directory (sess2analyse.curr_resdir).
%
% Input parameter:
%   S2A      struct with details of data to analyse (see getdata2analyse.m)

for snr = 1:length(s2a)
    
    try
        EEG = pop_loadset(fullfile(s2a(snr).folder,'EEG_2plot.set'));
    catch
        fprintf('No EEG\n');
        continue
        %         pause
    end
    
    EEG = pop_currentdensity(EEG, 'method','spline');
    
    try
        pop_saveset(EEG,fullfile(s2a(snr).folder,'EEG_2plot_CSD.set'))
    catch
        pause
    end
end
end

%--------------------------------------------------------------------------
function [EEG_ep1, EEG_ep2] = rej_badtrials(EEG_ep1,EEG_ep2,curr_resdir)

global  ALLEEG CURRENTSET EEG

% Reject bad trials/data
%-------------------------

EEG = EEG_ep2;

if EEG.trials~=1
    if exist([curr_resdir filesep 'rejected_epochs.mat'])==2; ifrej = 1; else; ifrej = 0; end;
% else
%     if exist([curr_resdir filesep 'events_with_boundaries.mat'])==2; ifrej = 1; else; ifrej = 0; end;
end


if ifrej==0
    pop_eegplot(EEG,1,1,1,[],'srate',EEG.srate,'spacing',75,...
        'eloc_file',EEG.chanlocs, 'winlength',30,'dispchans',32,'events',EEG.event,...
        'plottitle', 'Reject bad trials/data');
    
    
    fig1 = gcf;
    input('Select trials/data to reject, if ready press REJECT button, than any write any character to command window.\n');
end


% if EEG.trials~=1
    
    
    if ifrej==0
        allep = 1:size(EEG_ep1.epoch,2);
        rejected_eps = allep(~ismember(allep,[EEG.epoch.index]));
        save([curr_resdir filesep 'rejected_epochs.mat'],'rejected_eps');
    elseif ifrej==1
        load([curr_resdir filesep 'rejected_epochs.mat']);
        EEG = pop_rejepoch(EEG,rejected_eps);
    end
    
    EEG_ep1 = pop_rejepoch(EEG_ep1,rejected_eps);
    
    
    
% else
%     if ifrej==0
%         EEG_ep1.event = EEG.event;
%         events_with_boundaries = EEG.event;
%         save(fullfile(curr_resdir, 'events_with_boundaries.mat'),'events_with_boundaries')
%     elseif ifrej==1
%         load([curr_resdir filesep 'events_with_boundaries.mat'])
%         EEG_ep1.event = events_with_boundaries;
%         EEG.event = events_with_boundaries;
%     end
%     
% end
EEG_ep2 = EEG;
end


%--------------------------------------------------------------------------
function [EEG_ep1] = reject_bad_ica(EEG_ep1,EEG_ep2, curr_resdir)

% Reject bad trials manually
% Perform independent component analysis on high-pass filtered EEG data
% (EEG_ep2), removes manually selected ICA components from original (not
% high-pass filtered) data (EEG_ep1).
% curr_resdir: data folder to save new EEG data structure with ICA
% components (EEG_ICA.set) + save selected components (gcompreject.mat)


global ALLCOM ALLEEG CURRENTSET  EEG

% data for ICA analysis (it has to be assigned to EEG variable for the eeglab to properly execute ica related functions/GUIs)

%% Reject bad trials
if ~contains(curr_resdir,'LFP')
    
    [EEG_ep1, EEG_ep2] = rej_badtrials(EEG_ep1,EEG_ep2,curr_resdir);
    
end


%% ICA
if ~contains(curr_resdir,'LFP') && length(EEG_ep1.chanlocs)>1
    
    if EEG.trials~=1
        ica_setnm = [curr_resdir filesep 'EEG_ICA.set'];
        crej_nm= [curr_resdir filesep 'gcompreject.mat'];
    else
        ica_setnm = [curr_resdir filesep 'EEG_ICA_continu.set'];
        crej_nm= [curr_resdir filesep 'gcompreject_continu.mat'];
    end
    
    if exist(ica_setnm)==2; ifica = 1; else; ifica = 0; end;
    
    if ifica==0
        
        % get rank of data
        
        % curr_rank = rank(reshape(EEG_ep2.data,[size(EEG_ep2.data,1),size(EEG_ep2.data,2)*size(EEG_ep2.data,3)]));
        
        EEG = EEG_ep2;
        EEG_ICA = pop_runica(EEG, 'icatype', 'fastica');
        pop_saveset(EEG_ICA,ica_setnm);
        
    else
        EEG_ICA = pop_loadset(ica_setnm);
    end
    EEG = EEG_ICA; % it has to be assigned to EEG variable for the eeglab to properly execute ica related functions/GUIs
    
    
    % Label components to reject
    if exist(crej_nm)~=2
        try
            pop_eegplot( EEG, 0, 1, 1,[],'dispchans',20);
            EEG= pop_selectcomps(EEG, 1:20 );
        catch
            fprintf('ICA gone wild.\n')
            close(gcf); close(gcf);
            %continue
        end
        input('Select ICs to reject, if ready, press any key.\n');
    end
    %% Remove selected ICA components from original EEG
    EEG_ep1 = applyica(EEG,EEG_ep1,crej_nm);
    
end
end



%--------------------------------------------------------------------------
function EEG_ep1 = applyica(EEG,EEG_ep1,crej_nm)

if exist(crej_nm)~=2
    gcompreject = EEG.reject.gcompreject;
    save(crej_nm,'gcompreject');
else
    load(crej_nm);
end

% Apply ICA for "minimally" filtered data (dataset 1)
EEG_ep1.reject = EEG.reject; EEG_ep1.icawinv = EEG.icawinv;
EEG_ep1.icasphere = EEG.icasphere; EEG_ep1.icaweights = EEG.icaweights; EEG_ep1.icachansind = EEG.icachansind;
EEG_ep1.reject.gcompreject = gcompreject;

EEG_ep1 = pop_subcomp(EEG_ep1,[],1,0);
end


%--------------------------------------------------------------------------
function [EEG, Lnfilt] = rem_line_noise(EEG)
% Removes power line noise (50 Hz) from EEG data (eeglab structure)
% First tries CleanLine, if noise has not been removed sufficiently (has to
% be checked visually on the appeared PSD), notch filter is applied (45-55
% Hz notch filter). Label of applied filter is stored in Lnfilt variable.

% Is there any line noise?
[fig1, fig2] = check_eegdata(EEG);
inp0 = input('Linenoise? If no linenoise, press 0, otherwise any key.\n');
% inp0 = 1;
%
close(fig1); close(fig2);


if inp0==0
    Lnfilt = 'NoLN';
    
else
    % CLEANLINE FILTER
    
    EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:EEG.nbchan] ,...
        'computepower',1,'linefreqs',50,'newversion',0,...
        'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,...
        'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',4);
    Lnfilt = 'CLN';
    
    [fig1, fig2] = check_eegdata(EEG, strcat('CleanLine filtered data'));
    
    inp = input('Linenoise removed? If yes, press 1, otherwise any key\n');
    %     inp = 0;
    
    if inp~=1
        close(fig1,fig2)
        % notch filter to remove 50 Hz linenoise
        EEG = pop_eegfiltnew(EEG, 'locutoff',45,'hicutoff',55,'revfilt',1,'plotfreqz',1);
        Lnfilt = 'notch';
        [fig1, fig2] = check_eegdata(EEG,'');
    end
    
    close(fig1,fig2);
    
    
end
end

%--------------------------------------------------------------------------
function [fig1, fig2] = check_eegdata(EEG, figtitle)
% Plots power spectrum

narginchk(1,2)
if nargin<2
    figtitle = '';
end

% plot eeg time series

pop_eegplot(EEG,1,1,1,[],'srate',EEG.srate,'spacing',75,...
    'eloc_file',EEG.chanlocs, 'winlength',5,'dispchans',32,'events',EEG.event,...
    'plottitle', figtitle);
fig1 = gcf;

% plot power spectrum

fig2 = figure;
pop_spectopo(EEG,1,[1,size(EEG.data,2)],'EEG','freqrange',[1 100],'title',figtitle);
end


%--------------------------------------------------------------------------
function EEG_ep = substr_bas(EEG_ep,baseline_win)
% Subtract baseline from each channel and epoch

% EEG_ep: epoched EEG data (eeglab format)
% baseline_win: baseline window relative to event timetamps in seconds (ex: [-2 -1])
%-----------------------------------------------
bas_fr = abs(baseline_win(1,1)) + abs(baseline_win(1,2))*EEG_ep.srate;
EEG_ep.data = rmbase(EEG_ep.data,[],[1:bas_fr]); % dataset 1
end




%--------------------------------------------------------------------------
function fnm = find_filetag(allfiles,session,tag)


if (strcmp(session(end), 'l')||contains(session, 'left')) &&strcmp(tag, 'stimoff') %if left side - find the right marker file and behavior
    strmk = {'01' 'stim_off' 'stimoff' 'off'}; %stim off
elseif (strcmp(session(end), 'l')||contains(session, 'left'))&&strcmp(tag, 'stimon')
    strmk = {'03' 'stim_on' 'stimon' 'on'}; % stim on
elseif (strcmp(session(end), 'r')||contains(session, 'right'))&&strcmp(tag, 'stimoff') % right side
    strmk = {'02' 'stim_off' 'stimoff' 'off'}; %stim off
elseif (strcmp(session(end), 'r')||contains(session, 'right'))&&strcmp(tag, 'stimon') % right side
    strmk = {'04' 'stim_on' 'stimon' 'on'}; % stim on
end

for i = 1:length(allfiles)
    current_file = allfiles(i).name;
    if any(cellfun(@(x) contains(current_file,x), strmk))
        fnm = current_file;
    end
end
end