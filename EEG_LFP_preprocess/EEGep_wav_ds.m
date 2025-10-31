function [epoch_wavDS, epoch_wav,f] = EEGep_wav_ds(EEG_ep_choi,oc)
% EEGEP_WAV_DS Wavelet decomposition and downsampling of epoched data.
%   EEGep_wav_ds(EEG_ep_choi,oc) Performs wavelet decomposition on a single channel of 
%   epoched EEG/ LFP data. Epochs are concatenated and z-scored before wavelet decomp.
%   Wavelet coefficient matrix is downsampled. 
%   Input: 
%       EEG_EP_CHOI - eeglab data struct, contains epoched EEG/ LFP data
%       OC - double, index of channel to analyse
%
%   Output: 
%       EPOCH_WAVDS - downsampled 3D matrix (frequency x time x epochs)
%       EPOCH_WAV   - original 3D matrix (frequency x time x epochs)
%       F - double array, freq. components
%
% Johanna Petra Szabó, 10.2025
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


ds = floor(linspace(1,size(EEG_ep_choi.data,2),200));

% Concat. epoched data.
signconcat =  reshape(EEG_ep_choi.data(oc,:,:),1,[]);

% Zscore
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
end