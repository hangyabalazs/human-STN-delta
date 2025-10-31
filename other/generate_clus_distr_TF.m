function generate_clus_distr_TF(sess2analyse,epoch_win,baseline_win,varargin)
% GENERATE_CLUS_DISTR_TF Generates cluster distribution for cluster correction.
%   generate_clus_distr_TF(sess2analyse,epoch_win,baseline_win,...)
%   tests significant signal (EEG/ LFP/ WAVELET COHERENCE) change relative to baseline on multiple iterations.
%   Permutation test with FDR correction is applied on random data epochs selected from 
%   data of patients defined in SESS2ANALYSE. Significant clusters are selected.  
%   The cluster with the largest value (value defined by CLTYPE) is collected on each iteration for the predefined freq. bands, 
%   generating a distribution of significant clusters under the null for each freq. band separately, used
%   to identify appropriate threshold for cluster based correction (to
%   address multiple comparison).
%
% Required inputs:
%     SESS2ANALYSE      struct containing all necessary information (name of patient, side
%                       of experiment, tag of condition, session folder path) of
%                       session data that need to be analysed (see getdata2analyse)
%
%     EPOCH_WIN         1x2 vector, time window relative to event timestamp in sec, for data epoching (ex: [-2 2])
%
%     BASELINE_WIN       1x2 vector, time window relative to event timestamp in sec, for baseline correction
%
%
% Optional inputs (name-value pairs)
%   'stat_time'     1x2 vector, time window relative to event timestamp in sec, to perform statistics
% 
%   'freqs'         nx2 matrix, n = number of frequency bands to generate cluster distribution
% 
%   'evnr'          double, number of random data epochs
% 
%   'iter'          double, number of iterations to perform testing (= nr.of clusters within the distribution)
%
%   'evnr_stat'     double, number of epochs to include in stat. testing in
%                       one iteration, selected from EVNR nr of random epochs
%
%   'csd'           true | false, if true, CSD transformed EEG data is used
%                   (relevant only for postoperative EEG data)
% 
%   'bipol'         true | false, if true, F4-F3 bipolar derivation of EEG data is used
%                   (relevant only for postoperative EEG data)
% 
%   'choi'          cell array of channels of interest
%   
%   'CLtype'        string array, definition of cluster value: 
%                   'maxsum' - sum of t-values within the cluster
%                   'max'    - size of the cluster
%
%   'ALPHAS'        double or vector of alpha levels to use
%
%   'DATATYPE'      string array, type of data
%                   'potential' - EEG or LFP
%                   'coherence' - EEG-LFP wavelet coherence
%   'chanmean'      if true, channel averaged LFP data is used, if false
%                   LFP data is plotted channel-by-channel
%                   (relevant only for intraop LFP data)
%                   1 | 0 (default value: 1)
%
% See also: BOOSTAT_EEGLAB_J
%
% Johanna Petra Szabó, 10.2025
% Clinic for Neurosurgery and Neurointervention
% Semmelweis University, Budapest, Hungary
% szabo.johanna.petra@semmelweis.hu


prs = inputParser;
addRequired(prs,'sess2analyse',@isstruct);
addRequired(prs,'epoch_win',@isvector);
addRequired(prs,'baseline_win',@isvector);
addParameter(prs,'stat_time',[-1 1],@isvector);
addParameter(prs,'freqs',[.9 4.1; 3.9 7.1; 6.9 13.1; 12.9 30.1; 29.9 50.1; 49.9 80.1; 79.9 130.1],@ismatrix);
addParameter(prs,'evnr',100,@isnumeric);
addParameter(prs,'iter',200,@isnumeric);
addParameter(prs,'evnr_stat',32,@isnumeric);
addParameter(prs,'csd',true,@islogical);
addParameter(prs,'bipol',false,@islogical);
addParameter(prs,'choi',{},@iscell);
addParameter(prs,'CLtype','maxsum',@ischar);
addParameter(prs,'alphas',NaN,@(x) isnumeric(x)||isvector(x));
addParameter(prs,'datatype','potential',@ischar); % potential | coherence
addParameter(prs,'chanmean',false,@islogical);
parse(prs,sess2analyse,epoch_win,baseline_win,varargin{:});
p = prs.Results;

dbstop if error
global rootdir

if exist(fullfile(rootdir,'Cluster_distrib.mat'))~=2
    Cluster_distrib = struct;
else
    load(fullfile(rootdir,'Cluster_distrib.mat'))
end

rectime = p.sess2analyse(1).rectime;
rectype = p.sess2analyse(1).rectype;
p.rectype = rectype;

if ~strcmp(rectime,'postop')
    p.csd = false; 
end

if strcmp(rectype,'EEG')
    chanlocs = [];
    load(fullfile(rootdir,'postop_EEG_chanlocs.mat'));
    if p.bipol
        chanlabels = {'F4-F3'};
    else
        chanlabels = {chanlocs.labels};
    end
else
    chanlabels = arrayfun(@(x) ['Ch' num2str(x)] ,1:5,'UniformOutput',0);
end
if p.chanmean
    chnra = 1;
else
    chnra = length(chanlabels);
end
%%
epoch_wavDS = cell(length(sess2analyse),chnra );
for snr = 1:length(sess2analyse) % 19
    curr_resdir = sess2analyse(snr).folder;
    
    side = sess2analyse(snr).side;
    tag = sess2analyse(snr).tag;
    patnm = sess2analyse(snr).patient;
    
    if strcmp(p.datatype,'potential')
        [EEG_randE, ~] = randomEEG_epoch(curr_resdir, patnm, side, tag,epoch_win,p,[]);
        
        if isempty(EEG_randE); continue; end;
        %%
        
        
        chnr = min(length(chanlabels), size(EEG_randE.data,1));
        
        % Wavelet
        epochCh = cell(1,chnr);
        for oc = 1:chnr
            act_chan = chanlabels{oc};
            fprintf('%s ',act_chan);
            
            [epoch_wav, ~,f] = EEGep_wav_ds(EEG_randE,oc);
            
            epoch_p = abs(epoch_wav).^2;
            epochCh{oc} = epoch_p;
        end
       fieldnm = [rectime '_' rectype];
        
    elseif strcmp(p.datatype,'coherence')

       [epochCh, f] = random_wcoh(curr_resdir, patnm, side, tag,epoch_win,p) ;
       if isempty(epochCh); continue; end;
       
       fieldnm = 'EEG_LFP_wcoh';
    end
    
    if p.chanmean
        epoch_wavDS{snr,1} = mean(cat(4,epochCh{:}),4);
    else
        epoch_wavDS(snr,1:chnr) = epochCh;
    end
     
end

%%
formula = 'mean(arg1,3);';

times = epoch_win(1):1/50:epoch_win(2);
statlims = dsearchn(times',p.stat_time'); stattimes = times(statlims(1):statlims(2));
baslims = dsearchn(stattimes',p.baseline_win'); basinx = baslims(1):baslims(2);

exactp_c = cell(p.iter,chnra); maskersp = cell(p.iter,chnra);

for oc = 1:chnra
    if p.chanmean&&strcmp(rectype,'LFP')
        act_chan = 'chanmean';
    else
        act_chan = chanlabels{oc};
    end
    
    max_clustersize= nan(p.iter,size(p.freqs,1));
    for j = 1:p.iter % nr. of permutations
        randX = sort(randperm(p.evnr,p.evnr_stat));
        onech = epoch_wavDS(:,oc);
        nemp = find(~cellfun(@isempty,onech));
        onech_dat = cell(size(onech));
        onech_dat(nemp) = cellfun( @(x) x(:,:,randX) , onech(nemp),'UniformOutput',0) ;
        
        
        [bn, st, onech_datN, onech_avgs] = deal(cell(size(onech_dat)));
        bn(nemp) = cellfun(@(x) repmat(mean(x,2),[1 size(x,2) 1]), onech_dat(nemp),'UniformOutput',0);
        st(nemp) = cellfun(@(x) repmat(std(x,[],2),[1 size(x,2) 1]), onech_dat(nemp),'UniformOutput',0);
        onech_datN(nemp)  = arrayfun(@(x) (onech_dat{x}-bn{x})./st{x}, nemp,'UniformOutput',0);
        onech_avgs(nemp) = cellfun(@(x)  mean(x,3) ,onech_datN(nemp),'UniformOutput',0);
        dat = cat(3,onech_avgs{:});
        
        [exactp_ersp,maskersp{j,oc},alphafdr] = boostat_eeglab_J(dat(:,statlims(1):statlims(2),:),f,0.05,1000,false,'fdr',formula,basinx);
        p_c=  exactp_ersp;
        z_c = norminv(exactp_ersp);
        nsi = p_c>alphafdr;
        p_c(nsi) = 0;
        exactp_c{j,oc}= p_c;
        z_c(nsi)= 0;
        
        
        for k = 1:size(p.freqs,1)
            
            frix = f>=p.freqs(k,1)&f<p.freqs(k,2);
            act_islands = bwconncomp(p_c(frix,:));
            if numel(act_islands.PixelIdxList)>0
                
                % generating distribution according to largest clustersizes
                if strcmp(p.CLtype,'max')
                    act_clustsizes = cellfun(@length,act_islands.PixelIdxList);
                    max_clustersize(j,k) = max(act_clustsizes);
                elseif strcmp(p.CLtype,'maxsum')
                    % generating distribution according to largest sum of z-values included in clusters
                    
                    act_clustsizes = cellfun(@(x) sum(abs(z_c(x))),act_islands.PixelIdxList); % sum z values included in clusters
                    max_clustersize(j,k) = max(act_clustsizes); % find biggest sum
                end
            end
        end
        
    end
    
    Cluster_distrib.(fieldnm).(act_chan) = cell(size(p.freqs,1),3);
    for k = 1:size(p.freqs,1)
        Cluster_distrib.(fieldnm).(act_chan){k,1} = p.freqs(k,:);
        Cluster_distrib.(fieldnm).(act_chan){k,2} = max_clustersize(:,k);
        Cluster_distrib.(fieldnm).(act_chan){k,3} = prctile( max_clustersize(:,k) ,100-(100*0.05));
        
    end
    Cluster_distrib.Cell_column_labels = {'Freq_limits',[p.CLtype '_cluster_distribution'],'95th percentile'};
    save(fullfile(rootdir,'Cluster_distrib.mat'),'Cluster_distrib');
end
end


function [wcohE, ff] = random_wcoh(curr_resdir, patnm, side, tag,epoch_win,p)

[wcoh,f,ff,coi,chnames,phlag, wcohE]  = deal([]);
%% EEG


% Load EEG data
[EEG, randE] = randomEEG_epoch(curr_resdir, patnm, side, tag,epoch_win,p,[]);
if isempty(EEG); return; end;
eeg_chan = 'F4'; % ref: F3


%% LFP
% Load LFP data
curr_resdir_lfp =[curr_resdir '_LFP'];
[LFP,~] = randomEEG_epoch(curr_resdir_lfp, patnm, side, tag,epoch_win,p,randE);
if isempty(LFP); return; end;

chnames = {LFP.chanlocs.labels};


%%

% Prealloc.
[f,coi] = deal(cell(length(randE),1));

% Zscore data (across all epochs)
eegdats = permute(double(EEG.data),[ 2 3 1]);
eegdats_n = reshape( zscore( eegdats(:) ), size(eegdats));

wcohE= cell(1,length(chnames));
for ci = 1:length(chnames);
    
    
    lfpdats = permute(double(LFP.data(ci,:,:)),[ 2 3 1]);
    
    lfpdats_n = reshape( zscore( lfpdats(:) ), size(lfpdats));
    
    % Loop over events 
    
    [wcoh, phlag] = deal(cell(length(randE),1));
    for ei = 1:length(randE)
        eegdat = eegdats_n(:,ei); % one EEG epoch
        lfpdat = lfpdats_n(:,ei); % one LFP epoch
        
        % WAVELET COHERENCE
        [w_fall,wcross_fall,fall,coi{ei}] = wcoherence(eegdat,lfpdat,EEG.srate);
        
        %         % Phase Lag
        %         phL = angle(wcross_fall);
        
        % Get wav-coh within spec. freq. limits
        freqlim = [0.5 121];
        frinx = fall>=freqlim(1)&fall<=freqlim(2);
        wcoh{ei} = w_fall(frinx,:);
        %         phlag{ei} = phL(frinx,:);
        f{ei} = fall(frinx);
        
    end
    wcohE{ci} = cat(3,wcoh{:});
end
ff = f{1};
end



function [EEG_randE, randE] = randomEEG_epoch(curr_resdir, patnm, side, tag,epoch_win,p, randE)

    EEG_randE = [];
    fprintf('%s %s %s...\n',patnm, side, tag);
    try
        eegfnm = dir(fullfile(curr_resdir,['EEG_*_05HP.set']));
        if length(eegfnm); eegfnm = eegfnm(1); end
        EEG = pop_loadset(fullfile(curr_resdir,eegfnm.name));
        
    catch
        fprintf('NO EEG data %s %s %s...\n',patnm, side, tag);
        return;
    end
    
    if strcmp(p.rectype,'EEG')
        if p.csd
            
            EEG = pop_currentdensity(EEG, 'method','spline');
            
        elseif p.bipol
            EEG = pop_select(EEG,'channel',{'F4','F3'});
     
        end
    end
    
    % Generate timestamps for random events (skipping the first and last 3 secs)
    if isempty(randE)
        randE = sort(randperm([(EEG.pnts-6*EEG.srate)],p.evnr));
        randE = randE+ 3*EEG.srate;
    end
    
    EEG.event = struct;
    for k = 1:length(randE)
        EEG.event(k).latency =randE(k);
        EEG.event(k).type = 'Random';
        
    end
    
    % Epoching
    EEG_randE = pop_epoch(EEG, {}, epoch_win); % for analyses
end



