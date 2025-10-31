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
%
% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

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


c_min = [];
c_max = [];
for oc = 1:EEG_ep_choi.nbchan
    fprintf('Channel');
    act_chan = EEG_ep_choi.chanlocs(oc).labels;
    fprintf('%s ',act_chan);
    
    if exist([tf_dir filesep 'epoch_wav_ch' act_chan '_4.mat'])~=2
        
        
        [epoch_wavDS, epoch_wav,f] = EEGep_wav_ds(EEG_ep_choi,oc);
        
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