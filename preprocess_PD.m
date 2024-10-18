function preprocess_PD(sess2analyse,epoch_win,EventTypes_all,SubEventTypes_all,preproc,baseline_win)

%{
preprocess_PD(...) preprocesses EEG/LFP data for further analysis

Preprocessing steps:
1. Load raw data, import channel data, convert into EEGLAB dataset, save
             raw data is .set format.
2. Clean DBS artifacts (only if 'stimon' condition):
             DBSFILT GUI (by Guillaume Lio,2012, https://github.com/guillaumelio/DBSFILT)
             Applies frequency-domain Hampel filtering to clean DBS induced
             artifacts.
3. Import behavioral data into EEG dataset (according to TrialEvents.mat struct
             synchronized with EEG previously)
4. Resample: 250 Hz
5. Filter:   Lowpass filter with 100 Hz cutoff freq.
             CleanLine filter (by Tim Mullen, 2011, https://github.com/sccn/cleanline)
             Removes 50 Hz linenoise (uses sliding window to adaptively estimate sine wave amplitude to subtract)
             If not sufficiently effective (has to be approved manually),
             45 - 55 Hz notch filter is applied (pop_eegfiltnew).
             Two types of highpass filter (resulting in two datasets):
                - 0.5 Hz cutoff freq - for data used for further analysis
                - 2 Hz cutoff freq - for ICA
             Filtered datasets are saved in patient's result directory
             (curr_resdir in sess2analyse struct)
6. Data epoching (only if 'epochs' preprocessing is set!): data is epoched into trials, according to behavioral events
             and time window defined by epoch_win. Trial indeced are also saved
             separately in Evinxx.mat file.
7. Artifact rejection: - bad data portions/trials are removed manually from continuous/epoched EEG (pop_eegplot),
                            -- 'continu' preprocessing: boundaries (index & duration of rejected data) are inserted to EEG.event field
                            new EEG.event structure is saved in events_with_boundaries.mat file
                            -- 'epochs' preprocessing: indeces of rejected trials are saved into rejected_epochs.mat file
                       - (bad channels are removed and interpolated)
                       - ICA is performed following bad data/trial rejection
                       using 2Hz highpass filtered data
                       - Resulting components have to be reviewd visually.
                       Then components corresponding to artifacts related to
                       blinks/facial movements/bad channels are selected.
                       Selected components to reject are also saved separately
                       to 'gcompreject_continu.mat'
                       - Selected components are subtracted from the 0.5 Hz
                       filtered data.
9. Save final data structure: EEG_continu.set/EEG_2plot.set is used for further analysis.
                              Channel data also saved into a separate file:
                              [rectype '_EEG_chanlocs.mat']
       
Required inputs:
    sess2analyse - struct containing all necessary information (name of patient, side
                   of experiment, tag of condition, session folder path) of
                   sessions that need to be analysed (see getdata2analyse)
    epoch_win -  1x2 vector trial time window, relative to events in sec (ex: [-2 2])
    EventTypes_all - 1xN cell array of event labels ({'StimulusOn','Feedback','StopSignal','KeyPress1'};)
    SubEventTypes_all - Nx2 cell array of partition ("subevent") labels, each row
                    corresponds to an event label, columns to partitions
                    {'CueStim','StopStim';'Correct','Error';'FailedStopTrial','SuccesfulStopTrial';'CueResponse','StopResponse'};
    preproc - 'continu'| 'epochs' - preprocess and save eeg in coninuous/
                  epoched form

20-09-2022 Szabó Johanna Petra
Institute of Experimental Medicine, Budapest

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
        curr_resdir,EventTypes_all,SubEventTypes_all,rectype, rectime);
    
    
    if iste==0 % if is TrialEvents.mat (if synchronization was successful)
%         continue
    end
    
    if strcmp(preproc,'epochs')
        % Prepare epochs
        [EEG_filt1,EEG_filt2] =  prep_epochs(EEG_filt1, EEG_filt2, epoch_win);
        eeg_setnm = 'EEG_2plot.set';
        
        
    elseif strcmp(preproc,'continu')
        eeg_setnm = 'EEG_continu.set';
    end
    
    
    %     Reject bad data + apply ICA
    
    if exist(fullfile(curr_resdir,eeg_setnm))~=22
        [EEG] = reject_bad_ica(EEG_filt1,EEG_filt2, curr_resdir);
        
        
        if strcmp(preproc,'epochs')
            % Save event indeces
            save_evinxx(EEG,EventTypes_all,SubEventTypes_all,curr_resdir,true);
            
            
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
function [EEG_filt1 EEG_filt2 iste] = load_filter_data(currsess,patnm,side,condition,curr_resdir,EventTypes_all,SubEventTypes_all,rectype,rectime)



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
        EEG = behav_events(EEG,EventTypes_all,SubEventTypes_all,currsess,condition);
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
    [EEG Lnfilt] = rem_line_noise(EEG);
    
    
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


function reref2bipol(sess2analyse) 

% Re-reference data to bipolar montage (F3-F4)
% Saves re-referenced EEG data as EEG_2plot_bipol.set

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