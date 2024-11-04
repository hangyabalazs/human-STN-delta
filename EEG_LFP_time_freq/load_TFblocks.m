function [epoch_p f] = load_TFblocks(curr_resdir,act_chan,powphas,tf_dir,blocknr);
%LOAD_TFBLOCKS  Loads wavelet coefficients of epoched channel data
%   [epoch_p f] = LOAD_TFBLOCKS(curr_resdir,act_chan,powphas,tf_dir,blocknr)
%       Loads time-frequency matrices of complex wavelet coefficients derived from  
%       data of one channel (ACT_CHAN), saved in TF_DIR (a subdirectory of the patient's result directory, CURR_RESDIR).
%       Data is stored in BLOCKNR number of blocks (for faster load of data). 

% for each session and channel (1D: frequency, 2D: time, 3D: trials)

% Inputs parameters: 
%
%   CURR_RESDIR    path to results folder of patient's current session
%
%   ACT_CHAN       character array, label of channel of interest
%
%   POWPHAS        character array
%       'Pow'   outputs power coefficients
%       'Phas'  outputs phase coefficient       
%
%   TF_DIR         path to wavelet decomp results folder
%
%   BLOCKNR        double, nr of blocks used to save wavelet decomp results
%
% Output parameters:
%   EPOCH_P        time-frequency map (frequency x time x epochs), power or
%                  phase coefficients, depending on the value of POWPHAS
%
%   f              vector of frequency components
%
% See also: TIME_FREQ_PATIENTS

% Johanna Petra Szabó, Hangya Balázs, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu



narginchk(4,5)
if nargin<5
    blocknr = 4;
end

epoch_wav = [];
for bi = 1:blocknr
    % Load TF data
    try
        load(fullfile(tf_dir,['epoch_wav_ch' act_chan '_' num2str(bi) '.mat']));
        try
            load(fullfile(tf_dir,['epoch_f_ch' act_chan '.mat']));
        catch
            f = [];
        end
        
    
    ds = floor(linspace(1,1000,200));
    if size(epoch_wav_bl,2)>200
        epoch_wav_bl = epoch_wav_bl(:,ds,:); 
        save(fullfile(tf_dir,['epoch_wav_ch' act_chan '_' num2str(bi) '.mat']),'epoch_wav_bl')
    end
    
    catch
        try
            load(fullfile(curr_resdir,'Intraop_LFP','TFpows_blocks',['epoch_wav_ch' act_chan '_' num2str(bi) '.mat']));
            load(fullfile(curr_resdir,'Intraop_LFP','TFpows',['epoch_f_ch' act_chan '.mat']));
        catch
            fprintf('TFpow fail: %s.\n',curr_resdir);
            epoch_p = []; f = [];
            continue
        end
    end
    epoch_wav = cat(3,epoch_wav,epoch_wav_bl);
    epoch_wav_bl = [];
end


% Power coefficients
if strcmp(powphas, 'Pow')
    epoch_p = abs(epoch_wav).^2;
elseif strcmp(powphas, 'Phas')
    epoch_p = angle(epoch_wav);
end