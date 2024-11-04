function [wcoh,phlag,f,coi,chnames,TEinx_eeg2] = get_patient_wcoh(patnm,side,curr_resdir_eeg,event,subevent,downsamp);
% GET_PATIENT_WCOH  Wavelet coherence of intraop EEG and LFP in one patient
%   GET_PATIENT_WCOH(patnm,side,curr_resdir_eeg,event,subevent,downsamp) 
%       calculates mean squared wavelet coherence for each epoch around EVENT 
%       , partitioned according to SUBEVENT, of PATNM patient, SIDE side,
%       stored in CURR_RESDIR_EEG result folder. 
%
%   Input parameters:
%       PATNM               char. array of patient code
%       SIDE                char. array of tested side
%       CURR_RESDIR_EEG     result directory of EEG data
%       EVENT               char. array of event label
%       SUBEVENT            char. array of subevent label
%       DOWNSAMP            true | false, if true trial numbers are
%                           downsampled to match the subevent with the the least nr of trials 
%
%   Outout parameters:
%       WCOH                cell array of wavelet cross spectrum map
%       PHLAG               cell array of wavelet coherence phase lag maps
%       F                   cell array of frequency components
%       COI                 cell array with cone of influence in
%                           cycles/sample for the wavelet coherence
%       CHNAMES             cell array of channel labels
%       TEINX_EEG2          indeces of trials in 'TrialEvents*.mat'  file, 
%                           included in wavelet coherence, matched between 
%                           EEG and LFP data (epochs included in final
%                           preprocessed EEG/ LFP might slighly differ due to
%                           manual rejection of noisy trials) 
%
% See also: WCOHERENCE

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


[wcoh,f,coi,chnames,phlag]  = deal([]);
%% EEG

% Load EEG data
EEG = pop_loadset(fullfile(curr_resdir_eeg,'EEG_2plot.set'));
eeg_chan = 'F4'; % ref: F3

% Load event indeces for EEG data
load(fullfile(curr_resdir_eeg,'Evinxx.mat'));
TEinx_eeg = Evinxx.(event).(subevent).TE_index; % indeces of events correspoding to trials
evinx_eeg = Evinxx.(event).(subevent).epoch_index; % indeces of events correspoding to epochs in EEG data


%% LFP
% Load LFP data
curr_resdir_lfp =[curr_resdir_eeg '_LFP'];
LFP = pop_loadset(fullfile(curr_resdir_lfp,'EEG_2plot.set'));

% Load event indeces for LFP data
load(fullfile(curr_resdir_lfp,'Evinxx.mat'));
TEinx_lfp = Evinxx.(event).(subevent).TE_index; % indeces of events correspoding to trials
evinx_lfp = Evinxx.(event).(subevent).epoch_index; % indeces of events correspoding to epochs in LFP data

% if ~chanmean
    chnames = {LFP.chanlocs.labels};
% else
%     chnames = {'chanmean'};
% end

%%
% Indices of events existing in both EEG & LFP data (noisy epochs might have been rejected during preprocessing)
[~,ixa,ixb] = intersect(TEinx_eeg,TEinx_lfp);
evinx_eeg2 = evinx_eeg(ixa);
evinx_lfp2 = evinx_lfp(ixb);
TEinx_eeg2 = TEinx_eeg(ixa);
% Downsample trial nr to match contingency with minimum trial nr in order to compare contingencies (averaged wcoherence is affected by trial nr)
if downsamp && ismember(event,{'StimulusOn','StopSignal'})
    Len1 = length(Evinxx.(event).FailedStopTrial.epoch_index);
    Len2 = length(Evinxx.(event).SuccesfulStopTrial.epoch_index);
    Lenmin = min([Len1,Len2,length(evinx_eeg2)]);
    
    randev = sort(randperm(length(evinx_eeg2),Lenmin),'ascend');
    
    evinx_eeg2 = evinx_eeg2(randev);
    evinx_lfp2 = evinx_lfp2(randev);
end

if isempty(evinx_eeg2)
    fprintf('No trials for %s %s %s-%s\n',patnm,side,event,subevent)
    return
end
% Loop over LFP channels

% Prealloc.
[f,coi] = deal(cell(length(evinx_eeg2),1));
[wcoh, phlag] = deal(cell(length(evinx_eeg2),length(chnames)));

% Zscore data (across all epochs)
eegdats = permute(double(EEG.data),[ 2 3 1]);
eegdats_n = reshape( zscore( eegdats(:) ), size(eegdats));

for ci = 1:length(chnames);
    
%     if ~chanmean
        lfpdats = permute(double(LFP.data(ci,:,:)),[ 2 3 1]);
%     else
%         lfpdats = permute(double(mean(LFP.data,1)),[ 2 3 1]);
%     end
lfpdats_n = reshape( zscore( lfpdats(:) ), size(lfpdats));
    
    % Loop over events (existing in both kind of data)
    for ei = 1:length(evinx_eeg2)
        eegdat = eegdats_n(:,evinx_eeg2(ei)); % one EEG epoch
        lfpdat = lfpdats_n(:,evinx_lfp2(ei)); % one LFP epoch
        
        % WAVELET COHERENCE
        [w_fall,wcross_fall,fall,coi{ei}] = wcoherence(eegdat,lfpdat,EEG.srate);
        
        % Phase Lag
        phL = angle(wcross_fall);
        
        % Get wav-coh within spec. freq. limits
        freqlim = [0.5 80];
        frinx = fall>=freqlim(1)&fall<=freqlim(2);
        wcoh{ei,ci} = w_fall(frinx,:);
        phlag{ei,ci} = phL(frinx,:);
        f{ei} = fall(frinx);
        
        %              wcoherence(eegdat,lfpdat,EEG.srate,...
        %                 'FrequencyLimits',[0.5 100],'NumScalesToSmooth',20)
        %             pause(1); close(gcf)
    end
end




end