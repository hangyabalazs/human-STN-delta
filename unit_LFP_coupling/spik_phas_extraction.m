function spik_phas_extraction(sess2analyse,win,freqs, fr_names,PCdir,EventTypes,SubEventTypes,varargin)
%SPIK_PHAS_EXTRACTION  Phase values associated with spiking
% SPIK_PHAS_EXTRACTION(sess2analyse,win, freqs, fr_names,PCdir,EventTypes,SubEventTypes,...)
%   -Extracts phase values associated with spiking activity of units detected
%   in patients included in SESS2ANALYSE.
%   Phase values are extracted from time windows (PLOT_WIN) around a specific behavioral
%   events (EVENTTYPES) or subevents (SUBEVENTTYPES).
%   -Path to save results and plots is defined in PCDIR.
%
% Required inputs:
%     SESS2ANALYSE      struct containing all necessary information (name of patient, side
%                       of experiment, tag of condition, session folder path) of
%                       session data that need to be analysed (see getdata2analyse)
%     WIN          1x2 vector, time window relative to event timestamp in
%                       sec, to use for SPC calculation
%     PC_WIN            nr of smaller time windows (time window of WIN will
%                       be divided to PC_WIN nr of smaller windows to use for SPC calculation)
%     FREQS             Nx2 matrix, boundaries of frequency bands of
%                       interest (N is number of frequency bands to analyse, first column is
%                       the lower limit, second column is the lower limit)
%     FR_NAMES          cell array, labels of frequency bands of interest
%                       (nr of cells has to correspond to nr of rows in FREQS)
%     PCDIR             path to save results and figures
%
%     EVENTTYPES        1xN cell array of event labels, ex: {'StimulusOn','StopSignal','KeyPress1','Feedback'};
%
%     SUBEVENTTYPES     Nx2 cell array of partition ("subevent") labels, each row
%                       corresponds to an event label, columns to partitions
%                       {'CueStim','StopStim';'FailedStopTrial','SuccesfulStopTrial';'CueResponse','StopResponse';'Correct','Error';};
%
%    ISPLOT             1 | 0, if 1 figures are generated and saved
%
% Optional input (name-value pairs with default values):
%   'chanmean'      if true, channel averaged LFP data is used, if false
%                   LFP data is plotted channel-by-channel
%                   (relevant only for intraop LFP data)
%                   true | false (default value: false)
%   'subevs'        if true, PC is calculated for trials separated based on
%                   subevents deinfed in SUBEVENTTYPES
%                   true | false (default value: false)
%   'subregion'     character array or cell array, uses channel data
%                   derived from listed STN subregions
%                   (relevant only for intraop LFP data)
%          'all' | 'Motor' | 'Limbic' | 'Associative'  (default value: 'all')
%   'dominantfreq'  if true, dominant frequency within predefined frequency
%                   band limits (FREQS) are used for PC calculation, see FIND_DOMINANT_FREQ_BANDS
%                   true | false (default value: true)
%   'LFPdownsamp'   numeric value (N), downsample LFP signal by keeping every
%                   N-th sample starting with the first
%
% See also: FIND_DOMINANT_FREQ_BANDS

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


prs = inputParser;
addRequired(prs,'sess2analyse',@isstruct);
addRequired(prs,'win',@isvector);
addRequired(prs,'freqs',@ismatrix);
addRequired(prs,'fr_names',@iscell);
addRequired(prs,'PCdir',@isdir);
addRequired(prs,'EventTypes',@iscell);
addRequired(prs,'SubEventTypes',@iscell);
addParameter(prs,'chanmean',false,@islogical);
addParameter(prs,'subevs',false,@islogical);
addParameter(prs,'subregion','all',@ischar);
addParameter(prs,'dominantfreq',true,@islogical);
addParameter(prs,'LFPdownsamp',10,@isnumeric);
parse(prs,sess2analyse,win,freqs, fr_names,PCdir,EventTypes,SubEventTypes,varargin{:})
pr = prs.Results;




global rootdir
rectype = pr.sess2analyse(1).rectype;
rectime = pr.sess2analyse(1).rectime;
if ~strcmp(rectime,'intraop')
    error('Not intraop recordings, change rectime!');
end

frnr = length(pr.fr_names);




prev_curr_resdir = ''; prev_chan = '';

if strcmp(rectype,'LFP')&& pr.chanmean
    chlab = ['chanmean_' pr.subregion 'STN']; chtit = chlab;
elseif strcmp(rectype,'LFP')&& ~pr.chanmean
    chtit = 'by-channel';
elseif strcmp(rectype,'EEG')
    chlab = 'F4'; chtit = chlab;
end

if pr.subevs; subs = [1 2]; else; subs = 1; end;

if ~strcmpi(pr.subregion,'all') && strcmp(rectype,'LFP')
    load(fullfile(rootdir,'STN_loc','STN_loc.mat'));
    subreg_names = STN_loc.subreg_names;
    sinx = find(strcmp(subreg_names,pr.subregion));
end
nc = 1;


%%
for snr = 1:length(pr.sess2analyse)
    curr_resdir = pr.sess2analyse(snr).folder;
    side = pr.sess2analyse(snr).side;
    tag = pr.sess2analyse(snr).tag;
    patnm = pr.sess2analyse(snr).patient;
    currsess = pr.sess2analyse(snr).sessfolder;
    
    sL = strfind(currsess,filesep); sessnm = currsess(sL(end)+1:end);
    
    cellids = findcell('rat',patnm,'session',sessnm); cellnum = length(cellids);
    
    
    
    
    for iC = 1:cellnum
        
        act_cellid = cellids{iC};
        
        fprintf('%s...\n',act_cellid);
        
        cellidF =act_cellid;  cellidF(ismember(cellidF,'.')) = '_';
        
        [~,~,chan,unit] = cellid2tags(act_cellid);
        
        
        % Load epoched spike data
        SP = loadcb(act_cellid,'EVENTSPIKES');
        
        
        
        %% Load epoched EEG data
        if ~strcmp(curr_resdir,prev_curr_resdir) || (~strcmp(chan,prev_chan)&&~pr.chanmean&&strcmp(rectype,'LFP')) % units from the same patient (& channel) - spares loading again same LFP data
            
            
            
            % EEG/LFP signal
            if strcmp(rectype,'EEG')
                
                % Load preprocessed EEG
                EEG_ep = pop_loadset('filepath',curr_resdir,'filename','EEG_2plot.set');
                
                % Epoch indices
                load(fullfile(curr_resdir,'Evinxx.mat'));
                
                % Epoch timestamps
                
                orig_eeg_times = EEG_ep.times/1000;
                tinx = dsearchn(orig_eeg_times',pr.win');
                eeg_lims = orig_eeg_times(tinx);
                eeg_times = orig_eeg_times(tinx(1)+1:tinx(end));
                
                EEG_ep = pop_select(EEG_ep,'time',eeg_lims);
                
                
                
                
            elseif strcmp(rectype,'LFP')
                
                % Load raw LFP data
                EEG =  load_intraoplfp(currsess,patnm,side);
                
                
                % Load behav events
                EEG = behav_events(EEG,pr.EventTypes,pr.SubEventTypes,currsess,tag);
                
                % Downsample
                new_sr = round(EEG.srate/pr.LFPdownsamp);
                EEG = pop_resample(EEG,new_sr);
                
                
                % Split epochs
                [EEG_ep,~] =  prep_epochs(EEG, [], pr.win);
                
                
                % Make Evinxx struct (without saving)
                Evinxx =  save_evinxx(EEG_ep,EventTypes,SubEventTypes,curr_resdir,false);
                
                
                eeg_times = EEG_ep.times;
            end
        end
        
        prev_curr_resdir = curr_resdir;
        prev_chan = chan;
        
        % Select appropriate channel data (if not channel mean)
        if ~pr.chanmean && strcmp(rectype,'LFP')
            chlab = ['Ch' num2str(chan)];
            EEG_ep_ch = pop_select(EEG_ep,'channel',{chlab});
        else
            EEG_ep_ch = EEG_ep;
        end
        
        if isempty(EEG_ep_ch.data)
            fprintf('NO CHANNEL DATA FOR %s\',act_cellid)
            continue;
        end
        %% Find dominant frequency ranges within predefined frequency bands
        
        if pr.dominantfreq
            load(fullfile(rootdir,'DominantFreqs.mat'));
            
            
            newfreqs = nan(frnr,2);
            f_ind = cell(1,frnr);
            for frk = 1:frnr
                orig_frnm = fr_names{frk};
                
                try
                    newfreqs(frk,:) =  DominantFreqs.(orig_frnm).(patnm).(side).([rectime '_' rectype]).(tag).(chlab).new_fRange;
                    f_ind{frk} =DominantFreqs.(orig_frnm).(patnm).(side).([rectime '_' rectype]).(tag).(chlab).new_f_ind;
                catch
                    fprintf('No freq peak %s %s %s\n',patnm,chlab)
                    no_pk_cells{nc} = act_cellid;
                    nc = nc+1;
                    continue
                end
            end
        else
            newfreqs = pr.freqs;
        end
        
        %% Phase values derived from filtering + Hilbert transform
        hilb_phasdat = hilbEEG_epochs(EEG_ep_ch,newfreqs);
        %         hilb_phasdat = cellfun(@(x) mod(x,2*pi),hilb_phasdat0,'UniformOutput',0);
        
        
        %%
        
        for fri = 1:frnr
            
            if pr.dominantfreq && ~contains(pr.fr_names{fri},'dom')
                frnm = ['dom_' pr.fr_names{fri}];
            else
                frnm = [pr.fr_names{fri}];
            end
            
            
            
            for ei = 1:length(pr.EventTypes)
                event = pr.EventTypes{ei};
                
                
                % Find epochs in spike data
                evrow = ismember(SP.events(:,1),event);
                
                
                hilb_epochs = cell(1,length(subs));
                spike_win = cell(length(subs),1);
                epnr = nan(1,length(subs));
                for sei = subs
                    
                    if ~pr.subevs
                        evty = event;
                    else
                        evty = pr.SubEventTypes{ei,sei};
                    end
                    
                    % Find epochs in EEG data
                    
                    if ~strcmp(event,'StimulusOn') || (strcmp(event,'StimulusOn')&& ~pr.subevs)
                        evinx = Evinxx.(event).(evty).epoch_index;
                        TE_index = Evinxx.(event).(evty).TE_index;
                        
                        
                    elseif strcmp(event,'StimulusOn')&& pr.subevs
                        [TE_index,evinx] = StimOn_stoppart_evinx(Evinxx,event,evty);
                    end
                    
                    
                    
                    % Get epochs of interest
                    hilb_epochs{sei} = hilb_phasdat{fri}(1,:,evinx);
                    
                    spike_epochs = SP.event_stimes{evrow,1}(TE_index);
                    epnr(sei)= length(TE_index);
                    
                    
                    % Spike timestamps according to win
                    spike_win{sei,1} = cellfun(@(x) x(x>=win(1)& x<=win(2)), spike_epochs, 'UniformOutput',0);
                    
                    
                    
                end; clear sei;
                
                
                
                
                
                for sei = subs
                    
                    if ~pr.subevs
                        evty = event;
                    else
                        evty = pr.SubEventTypes{ei,sei};
                    end
                    
                    
                    resdir = fullfile(pr.PCdir,chtit,event,evty);
                    if ~isdir(resdir); mkdir(resdir); end;
                    
                    
                    
                    
                    
                    %% Indices of spikes
                    
                    
                    % sp_hilb = cellfun(@(x) interp_phasvals(x,eeg_times,hilb_epochs{sei}),spike_win(sei,:),'UniformOutput',0);
                    sp_hilb = cell(1,epnr(sei));
                    
                    for j = 1:epnr(sei)
                        spw = spike_win{sei,:}{j};
                        
                        if ~isempty(spw)
                            sp_inx = dsearchn(eeg_times',spw); % spike index in eeg/lfp signal
                            sp_hilb{j} =  hilb_epochs{sei}(1,sp_inx,j); % associated phase value (within one epoch)
                        end
                        
                    end
                    
                    sp_ts = cellfun(@(x) cat(1,x{:}), spike_win(sei,:),'UniformOutput',0);
                    
                    %%
                    
                    
                    
                    % Complete (instead of overwriting) results mat
                    if exist(fullfile(resdir,'PC_phase_values_allcells.mat'))==2
                        load(fullfile(resdir,'PC_phase_values_allcells.mat'))
                    else
                        PC_phase_values_allcells = struct;
                    end
                    
                    
                    
                    PC_phase_values_allcells.(cellidF).(frnm).Hilb_phase_values =  [sp_hilb{:}];
                    
                    PC_phase_values_allcells.(cellidF).(frnm).spike_timestamps = sp_ts{:};
                    
                    
                    
                    
                    % save matrix
                    save(fullfile(resdir,'PC_phase_values_allcells.mat'),'PC_phase_values_allcells');
                end
            end
        end
    end
end
end





%---------------------------------------------------------------------------
function hilb_phasdat = hilbEEG_epochs(EEG_ep,freqs)

frnr = size(freqs,1);

hilb_phasdat = cell(1,frnr);

for fr = 1:frnr
    
    % bandpass filter EEG/LFP
    EEGfilt = pop_eegfiltnew(EEG_ep,freqs(fr,1),freqs(fr,2));
    
    if size(EEG_ep.data,1)>1
        sign = double(mean(EEGfilt.data,1));
    else
        sign = double(EEGfilt.data);
    end
    
    sign_re = reshape(sign,[size(sign,1), size(sign,2)*size(sign,3)]);
    
    hilbsign = hilbert(sign_re);
    
    hilb_phasdat{fr} = reshape(angle(hilbsign),[size(sign,1), size(sign,2), size(sign,3)]);
end
end