function find_dominant_freq_bands(sess2analyse,fr_names,freq_bands,chanmean,subregion)
%FIND_DOMINANT_FREQ_BANDS   Selects dominant freuency within band limits
%   FIND_DOMINANT_FREQ_BANDS(sess2analyse,fr_names,freq_bands,chanmean,subregion)
%       - loads full-epoch normalized, epoched time-frequency data
%       - selects frequency components within band limits
%       - averages across all timepoints and all epochs
%       - selects frequency components with the highest power value (=dominant frequency)
%       - defines a frequency range around the dominant frequency (dom. freq.): dom. freq. +- dom. freq./6 Hz
%       - stores results of all patients in a struct, saved as 'DominantFreqs.mat' in rootdir
%
% See also: TIME_FREQ_PATIENTS, LOAD_TFBLOCKS

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

global rootdir

if exist(fullfile(rootdir,'DominantFreqs.mat'))~=2
    DominantFreqs = struct;
else
    load(fullfile(rootdir,'DominantFreqs.mat'));
end


rectime = sess2analyse(1).rectime;
rectype = sess2analyse(1).rectype;

if ~strcmp(subregion,'all') && strcmp(rectype,'LFP')
    load(fullfile(rootdir,'STN_loc','STN_loc.mat'));
    subreg_names = STN_loc.subreg_names;
    sinx = find(strcmp(subreg_names,subregion));
end
frnr = length(fr_names);


if strcmp(rectype,'EEG') && strcmp(rectime,'postop')
    load(fullfile(rootdir,'postop_EEG_chanlocs.mat')); choi = {chanlocs.labels}; clear chanlocs;
elseif strcmp(rectype,'EEG') && strcmp(rectime,'intraop')
    choi = {'F4'};
elseif strcmp(rectype,'LFP') && chanmean==1
    choi = {['chanmean_' subregion 'STN']};
elseif strcmp(rectype,'LFP') && chanmean~=1
    choi = {'Ch1','Ch2','Ch3','Ch4','Ch5'};
end



%%
for snr = 1:length(sess2analyse)
    curr_resdir = sess2analyse(snr).folder;
    side = sess2analyse(snr).side;
    tag = sess2analyse(snr).tag;
    patnm = sess2analyse(snr).patient;
    
    fprintf('%s...\n',curr_resdir)
    
    tf_dir = fullfile(curr_resdir,'TFpows_blocks');
    epoch_pow = [];
    
    for ci = 1:length(choi)
        act_chan = choi{ci};
        epoch_pow = [];
        % TF data
        if ~chanmean && strcmp(rectype,'LFP')
            
            [epoch_pow f] = load_TFblocks(curr_resdir,act_chan,'Pow',tf_dir,4);
            
        elseif chanmean==1 && strcmp(rectype,'LFP')
            tfs = dir(fullfile(curr_resdir,'TFpows_blocks','epoch_wav*_1.mat'));
            
            dat_4D_pow = [];
            for ti =1:length(tfs)
                anm = tfs(ti).name;
                lstr = strfind(anm,'_');
                ach = anm(strfind(anm,'Ch'):lstr(end)-1);
                
                if ~strcmp(subregion,'all')
                    if ~ismember(STN_loc.Patients.(patnm).(side).Channels.(ach).Subregion(sinx),[1 2]);
                        fprintf('Not %s: %s %s %s\n',subregion,patnm,side,ach)
                        continue
                    end
                end
                
                [e_pow f] = load_TFblocks(curr_resdir,ach,'Pow',tf_dir,4);
                dat_4D_pow =  cat(4,dat_4D_pow,e_pow);
            end
            if isempty(dat_4D_pow)
                fprintf('Not %s ATALL: %s %s \n',subregion,patnm,side)
                continue
            end
            epoch_pow = nanmean(dat_4D_pow,4);
            
        elseif strcmp(rectype,'EEG')
            [epoch_pow f] = load_TFblocks(curr_resdir,act_chan,'Pow',tf_dir,4);
            
        end
        
        
        %% Find dominant frequency ranges within predefined frequency bands
        
        
        % Full-epoch length single-trial normalization
        M = median(epoch_pow,2);
        Mrep = repmat(M,[1 size(epoch_pow,2)]);
        evpowTRN = epoch_pow./Mrep;
        
        for frk = 1:frnr
            
            orig_frnm = fr_names{frk};
            orig_freqs = freq_bands(frk,:);
            
            pre_f_inx = find(f>freq_bands(frk,1)&f<freq_bands(frk,2));
            
            if ~isempty(evpowTRN)
                psd = mean(evpowTRN(pre_f_inx,:,:),[2 3]);
                
                [pks,locs] = findpeaks(psd);
                
                
                if ~isempty(locs)
                    [~,maxipk] = max(pks);
                    maxilc = locs(maxipk);
                    dom_f_ind = pre_f_inx(maxilc);
                    dom_freq = f(dom_f_ind);
                    dom_frange = [dom_freq-(dom_freq/6) dom_freq+(dom_freq/6)];
                    dom_fra_ind = intersect(find(f>=dom_frange(1)),find(f<=dom_frange(2)));
                else
                    %                     fprintf('No freq peak %s %s\n',patnm,act_chan)
                    %                     DominantFreqs.(orig_frnm).(patnm).(side).([rectime '_' rectype]).(tag).(act_chan) = [];
                    %                     continue
                    dom_freq = median(f(pre_f_inx));
                    dom_frange = [dom_freq-(dom_freq/6) dom_freq+(dom_freq/6)];
                    dom_fra_ind = intersect(find(f>=dom_frange(1)),find(f<=dom_frange(2)));
                end
            else
                DominantFreqs.(orig_frnm).(patnm).(side).([rectime '_' rectype]).(tag).(act_chan) = [];
                continue;
                
            end
            
            DominantFreqs.(orig_frnm).(patnm).(side).([rectime '_' rectype]).(tag).(act_chan).dom_freq = dom_freq;
            DominantFreqs.(orig_frnm).(patnm).(side).([rectime '_' rectype]).(tag).(act_chan).new_fRange = dom_frange;
            DominantFreqs.(orig_frnm).(patnm).(side).([rectime '_' rectype]).(tag).(act_chan).orig_fRange = orig_freqs;
            
            DominantFreqs.(orig_frnm).(patnm).(side).([rectime '_' rectype]).(tag).(act_chan).new_f_ind = dom_fra_ind;
            DominantFreqs.(orig_frnm).(patnm).(side).([rectime '_' rectype]).(tag).(act_chan).orig_f_ind = pre_f_inx;
        end
    end
end

save(fullfile(rootdir,'DominantFreqs.mat'),'DominantFreqs')