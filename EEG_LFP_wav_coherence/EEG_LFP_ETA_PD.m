function EEG_LFP_ETA_PD(EventTypes)
%EEG_LFP_ETA_PD     Intraop EEG and LFP event-triggered averages (ETA) plotted aligned
%   EEG_LFP_ETA_PD  draws ETAs of EEG and LFP patient-by-patient and
%   averaged across patients.
%
%   Input parameter:
%       EVENTTYPES   cell array of event labels
%

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


global filesdir

% Parameters
plot_win = [-1 1];
baseline_win = [-1 -.5];
freqs = [1 4];

sess2analyse = getdata2analyse(filesdir, 'rectype','EEG',...
    'rectime','intraop','patients', 'allpatients', 'side','left', 'condition','stimoff');

% EEG + LFP aligned - patient-by-patient (plot + save struct)

for ei = 1:length(EventTypes)
    event = EventTypes{ei};
    
    eeg_lfp_suprapus(sess2analyse,event,{},freqs,plot_win,true,baseline_win,[],1)
end


% EEG + LFP aligned - averaged across patients (plot + save struct)



for ei = 1:length(EventTypes)
    event = EventTypes{ei};
    
    eeg_lfp_AVG_suprapus(event,'freqs',freqs,'plot_win',plot_win,'baseline_win',baseline_win);
    
end




end





%--------------------------------------------------------------------------
function eeg_lfp_suprapus(sess2analyse,event,subevs,freqs,plot_win,isfig,baseline_win,movavg_dp,chanmean)

global figdir_pd


srate = 250;
time = -2 :1/srate: 2; time = time(1:end-1);
pltime = plot_win(1):1/srate:plot_win(2);
timinx = dsearchn(time',pltime');



% Result directory + load result structs if they already exist
savedir = fullfile(figdir_pd,'introp_EEG_LFP','ETA');
if ~isempty(baseline_win)
    eeg_structname = fullfile(savedir,['EEG_ETA_freq' num2str(freqs) '_bas' num2str(baseline_win) '.mat']);
    lfp_structname = fullfile(savedir,['LFP_ETA_freq' num2str(freqs) '_bas' num2str(baseline_win) '.mat']);

else
    eeg_structname = fullfile(savedir,['EEG_ETA_freq' num2str(freqs) '.mat']);
    lfp_structname = fullfile(savedir,['LFP_ETA_freq' num2str(freqs) '.mat']);
end

if exist(eeg_structname)==2
    load(eeg_structname); else; EEG_ETA = struct;
end

if exist(lfp_structname)==2
    load(lfp_structname); else; LFP_ETA = struct;
end


if isempty(subevs); subnrs = 3; else; subnrs = 1:2; end;


% PATIENT LOOP
for k = 1:length(sess2analyse)
    
    patnm = sess2analyse(k).patient;
    side = sess2analyse(k).side;
    fprintf('%s %s ...\n',patnm,side);
    
    %% EEG
    
    curr_resdir_eeg = sess2analyse(k).folder;
    % Load EEG data
    EEG = pop_loadset(fullfile(curr_resdir_eeg,'EEG_2plot.set'));
    eeg_chan = 'F4'; % ref: F3
    
    
    %% LFP
    % Load LFP data
    curr_resdir_lfp =[curr_resdir_eeg '_LFP'];
    LFP = pop_loadset(fullfile(curr_resdir_lfp,'EEG_2plot.set'));
    
    
    %% Loop over subevents (if subevent are given)
    for ei = subnrs
        if ei==3;
            subevent = event; evdir = event;
        else;
            subevent = subevs{ei}; evdir = [event '_' subevent];
        end;
        
        % Result sub-directory
        if isempty(freqs)
            figdir = fullfile(figdir_pd,'introp_EEG_LFP','ETA',evdir);
        else
            frnm = ['F' num2str(freqs)];
            figdir = fullfile(figdir_pd,'introp_EEG_LFP','ETA',evdir,frnm);
        end
        if ~isfolder(figdir); mkdir(figdir); end;
        
        
        % Load event indeces for EEG data
        load(fullfile(curr_resdir_eeg,'Evinxx.mat'));
        TEinx_eeg = Evinxx.(event).(subevent).TE_index; % indeces of events correspoding to behavioral trials
        evinx_eeg = Evinxx.(event).(subevent).epoch_index; % indeces of events correspoding to epochs in EEG data
        
        
        % Load event indeces for LFP data
        load(fullfile(curr_resdir_lfp,'Evinxx.mat'));
        TEinx_lfp = Evinxx.(event).(subevent).TE_index; % indeces of events correspoding to behavioral trials
        evinx_lfp = Evinxx.(event).(subevent).epoch_index; % indeces of events correspoding to epochs in LFP data
        
        if ~chanmean
        chnames = {LFP.chanlocs.labels};
        else
            chnames = {'chanmean'};
        end
        
        %
        % Indices of events existing in both EEG & LFP data (noisy epochs might have been rejected during preprocessing)
        [~,ixa,ixb] = intersect(TEinx_eeg,TEinx_lfp);
        evinx_eeg2 = evinx_eeg(ixa); 
        evinx_lfp2 = evinx_lfp(ixb);
        eeg_chan = 'F4';
        
        if isempty(evinx_eeg2);
            fprintf('No %s trials for %s %s\n',subevent,patnm,side);
            continue;
        end
        
        
        % Loop over LFP channels
        for ci = 1:length(chnames);
            lfp_chan = chnames{ci};
            
            eegdats = permute(double(EEG.data(1,:,evinx_eeg2)),[2 3 1]); % EEG epochs
            
            if ~chanmean
                lfpdats = permute(double(LFP.data(ci,:,evinx_lfp2)),[2 3 1]); % LFP epochs
            else
                
                lfpdats = permute(double(nanmean( LFP.data(:,:,evinx_lfp2) ,1)),[2 3 1]); % LFP epochs
            end
            
            % Filter
            if ~isempty(freqs)
                eegdats0 = eegdats;
                eegdats = filtdat(eegdats0,freqs,srate);
                lfpdats0 = lfpdats;
                lfpdats = filtdat(lfpdats0,freqs,srate);
            end
            
            if ~isequal(time,pltime); isfig2 = false; else; isfig2 = true; end;
            
            % Save ETA in original time window
            [eegdats_avg, eegdats_std, lfpdats_avg, lfpdats_std] = ...
                eeg_lfp_etafig(time,eegdats,lfpdats,event,...
                patnm,side,lfp_chan,figdir,isfig2,baseline_win,[]);
            
            EEG_ETA.(patnm).(side).(eeg_chan).(event).(subevent).AVG = eegdats_avg;
            EEG_ETA.(patnm).(side).(eeg_chan).(event).(subevent).STD = eegdats_std;
            
            LFP_ETA.(patnm).(side).(lfp_chan).(event).(subevent).AVG = lfpdats_avg;
            LFP_ETA.(patnm).(side).(lfp_chan).(event).(subevent).STD = lfpdats_std;
            
            % Plot ETA in diff time window (with baseline norm. if applicable)
            if ~isfig2 && isfig
                eeg_lfp_etafig(pltime,eegdats(timinx,:),lfpdats(timinx,:),event,...
                    patnm,side,lfp_chan,figdir,true,baseline_win,movavg_dp);
            end
            
        end
    end
end

save(eeg_structname,'EEG_ETA')
save(lfp_structname,'LFP_ETA')
end





%--------------------------------------------------------------------------
function  eeg_lfp_AVG_suprapus(event,varargin)

prs = inputParser;
addRequired(prs,'event',@ischar)
addParameter(prs,'subevs',{},@(x) ischar(x)|iscell(x)) % {} | '' |
addParameter(prs,'patients','all',@(x) ischar(x)|iscell(x)) % 'all' | {'pd01','pd03',...}
addParameter(prs,'region',{},@iscell) % 'all' | {'Motor'} | {'Associative'} | {'Limbic'}
addParameter(prs,'close2centr','all',@ischar) % 'all' | 'Motor' | 'Associative' | 'Limbic'
addParameter(prs,'maxpower','all',@ischar)  % 'all' | 'high_delta' | 'beta' | 'fast_gamma'
addParameter(prs,'chanmean',1,@isnumeric)
addParameter(prs,'side','left',@(x) ischar(x)|iscell(x)) % 'both' | {'left','right'}
addParameter(prs,'iscell',[],@(x) islogical(x)| isnumeric(x)) % [] | true | false
addParameter(prs,'freqs',[],@(x) isvector(x)| isnumeric(x));
addParameter(prs,'isfig',true,@islogical);
addParameter(prs,'plot_win',[-2 2],@isvector);
addParameter(prs,'baseline_win',[],@(x) isvector(x)| isnumeric(x));
addParameter(prs,'movavg_dp',[],@(x) isvector(x)| isnumeric(x));

parse(prs,event,varargin{:})
g = prs.Results;

global rootdir figdir_pd


srate = 250;
time = -2 :1/srate: 2; time = time(1:end-1);
pltime = g.plot_win(1):1/srate:g.plot_win(2);
timinx = dsearchn(time',pltime');

savedir = fullfile(figdir_pd,'introp_EEG_LFP','ETA');

% Load presaved ETA for each patient
eeg_structname = fullfile(savedir,['EEG_ETA_freq' num2str(g.freqs) '_bas' num2str(g.baseline_win) '.mat']);
lfp_structname = fullfile(savedir,['LFP_ETA_freq' num2str(g.freqs) '_bas' num2str(g.baseline_win) '.mat']);
 
load(eeg_structname)
load(lfp_structname)


% Patient names
eeg_chan = 'F4';
if strcmp(g.patients,'all')
    pats = fieldnames(EEG_ETA);
else
    pats = g.patients;
end


% Patient sides
if strcmp(g.side,'both')
    sids = {'left','right'};
else
    if ~iscell(g.side); sids = {g.side}; else; sids = g.side; end;
end

% Event/ subevent iterations
if isempty(g.subevs); subnrs = 3; else; subnrs = 1:2; end;


for ei = subnrs
    if ei==3;
        subevent = event; evdir = event;
    else;
        subevent = g.subevs{ei}; evdir = [event '_' subevent];
    end;
    
    
    % FIGURE DIRECTORY to save
    if isempty(g.freqs)
        figdir = fullfile(figdir_pd,'introp_EEG_LFP','ETA',evdir);
    else
        frnm = ['F' num2str(g.freqs) '_BAS' num2str(g.baseline_win)];
        figdir = fullfile(figdir_pd,'introp_EEG_LFP','ETA',evdir,frnm);
    end
    if ~isfolder(figdir); mkdir(figdir); end;
    
    
    
    
    k= 1;
    [eegdats,lfpdats] = deal(cell(1,length(pats)*length(sids)));
    
    for pp = 1:length(pats)
        patnm = pats{pp};
%         if strcmp(patnm,'pd05'); continue; end;
        
        for ss = 1:length(sids)
            sidenm = sids{ss};
            
            %% Select channels to average
            if ~g.chanmean
                inchans = {};
                subreg_tit = '';
                % Get channels in predefined subregion
%                 if ~isempty(g.region{1});
%                     inchans0 = cell(1,length(g.region));
%                     for rr = 1:length(g.region)
%                         try
%                             inchans0{rr} = get_chan_subregion(patnm,sidenm,g.region{rr});
%                         catch
%                             disp('');
%                         end
%                     end
%                     inchans = cat(2,inchans0{:});
%                     subreg_tit = [g.region{:}];
%                 end
                
                
                
                % Get channel closest to centroid of one STN subreg.
                if  ~strcmp(g.close2centr,'all')
                    load(fullfile(rootdir,'Channels_close2centroids.mat'));
                    inchans = centr_tab{[patnm '_' sidenm(1)],g.close2centr};
                    subreg_tit = ['Close to ' g.close2centr ' centr.'];
                    
                end
                
                % Get channel with maximum power in spec. freq band
                if  ~strcmp(g.maxpower,'all')
                    load(fullfile(rootdir,'maxFREQ_channels.mat'));
                    inchans =  maxFREQ_channels{[patnm '_' sidenm(1)],g.maxpower};
                    subreg_tit = ['Max ' g.maxpower ' power'];
                    
                end
                
                % Get channels with cells
                if ~isempty(g.iscell)
                    cids = findcell('rat',patnm); % cellids for patient
                    cellids = cids(contains(cids,sidenm(1))); % cellids on the current side of patient
                    
                    cell_chan = nan(1,length(cellids));
                    for cc = 1:length(cellids)
                        [~,~,cell_chan(cc),~] = cellid2tags(cellids{cc}); % channels on which cells were detected
                    end
                    inchans0 = arrayfun(@(x) ['Ch' num2str(x)],unique(cell_chan),'UniformOutput',0);
                    
                    if ~g.iscell
                        allchans = fieldnames(LFP_ETA.(patnm).(sidenm));
                        
                        inchans00 = inchans0;
                        inchans0 = allchans( ~ismember(allchans,inchans00) );
                    end
                    
                    if isempty(inchans); inchans = inchans0; else; inchans = intersect(inchans,inchans0); end;
                    
                    cell_tit = ['cells-' char(string(g.iscell))];
                else
                    cell_tit = '';
                end
                
%                 % If no condition, get all channels
%                 if isempty(g.region{1}) && isempty(g.iscell) && strcmp(g.close2centr,'all') && strcmp(g.maxpower,'all')
%                     inchans = fieldnames(LFP_ETA.(patnm).(sidenm));
%                 end
                
            else
                inchans = {'chanmean'};
            end
            
            chnr = length(inchans);
            fprintf('%d nr of channels for %s, %s\n',chnr,patnm,sidenm);
            
            
            lfpavgs = cell(1,chnr);
            for ci = 1:chnr
                lfpchan = inchans{ci};
                if isfield(LFP_ETA.(patnm).(sidenm),lfpchan)
                    if isfield(LFP_ETA.(patnm).(sidenm).(lfpchan).(event),subevent)
                        lfpavgs{ci} = LFP_ETA.(patnm).(sidenm).(lfpchan).(event).(subevent).AVG;
                    else
                        continue;
                    end
                else
                    fprintf('No %s in %s %s\n',lfpchan,patnm,sidenm);
                end
            end
            if isfield(EEG_ETA.(patnm).(sidenm).(eeg_chan).(event),subevent)
                eegdats{k} = EEG_ETA.(patnm).(sidenm).(eeg_chan).(event).(subevent).AVG;
            else
                continue;
            end
            if ~isempty(cat(2,lfpavgs{:}))
                lfpdats{k} = mean(cat(2,lfpavgs{:}),2);
                k = k+1;
            end
            
        end
    end
    
    lfps = cat(2,lfpdats{:});
    eegs = cat(2,eegdats{:});
    
    if ~g.chanmean
        titL =  [subreg_tit ' ' cell_tit];
    else
        titL =  'chanmean';
    end
    
    eeg_lfp_etafig(time(timinx),eegs(timinx,:),lfps(timinx,:),event,'PatAVG',[g.side 'side'],...
       titL,figdir,true,[],g.movavg_dp);
    
end
end



%--------------------------------------------------------------------------
function [eegdats_avg, eegdats_std, lfpdats_avg, lfpdats_std] = ...
    eeg_lfp_etafig(time,eegdats,lfpdats,event,patnm,side,lfp_chan,figdir,isfig,baseline_win,movavg_dp)

Len = length(time);

% Apply baseline normalization
if ~isempty(baseline_win)
    baslim_inx = dsearchn(time',baseline_win'); basinx = baslim_inx(1):baslim_inx(2);
    
    % EEG basline normalized
    repbas_mean = repmat( mean(eegdats(basinx,:),1) ,[Len 1]); 
    repbas_std =  repmat( std(eegdats(basinx,:),[],1) ,[Len 1]);
    
    eegdats = ( eegdats - repbas_mean ) ./ repbas_std;
    
    % LFP basline normalized
    repbas_mean = repmat( mean(lfpdats(basinx,:),1) ,[Len 1]);
    repbas_std =  repmat( std(lfpdats(basinx,:),[],1) ,[Len 1]);
    
    lfpdats = ( lfpdats - repbas_mean ) ./ repbas_std;
    
end

% EEG avg & std
eegdats_n = reshape( zscore( eegdats(:) ), size(eegdats));
eegdats_std = std(eegdats_n, [], 2);
eegdats_std = eegdats_std / sqrt(size(eegdats_n, 2));
eegdats_avg = mean(eegdats_n,2);

% LFP avg & std
lfpdats_n = reshape( zscore( lfpdats(:) ), size(lfpdats));
lfpdats_std = std(lfpdats_n, [], 2);
lfpdats_std = lfpdats_std / sqrt(size(lfpdats_n, 2));
lfpdats_avg = mean(lfpdats_n,2);

% Moving average
if ~isempty(movavg_dp)
    eegdats_avg = movavg(eegdats_avg,'simple',movavg_dp);
    eegdats_std = movavg(eegdats_std,'simple',movavg_dp);
    lfpdats_avg = movavg(lfpdats_avg,'simple',movavg_dp);
    lfpdats_std = movavg(lfpdats_std,'simple',movavg_dp);
end

k = size(lfpdats,2);

if isfig
    fig = figure;
    
    h1 = errorshade(time,eegdats_avg,eegdats_std,'LineColor',[0 0 0.6],'ShadeColor',[0 0 0.6],'LineWidth',1.5);
    hold on;
    h2 = errorshade(time,lfpdats_avg,lfpdats_std,'LineColor',[0 0.6 0.6],'ShadeColor',[0 0.6 0.6],'LineWidth',1.5);
    xlim([time(1) time(end)])
    legend([h1(1) h2(1)],{'EEG','LFP'});
    xlabel(['Time relative to ' event ' (s)']);
    title({[event ' triggered average'],[patnm ' ' side '(n=' num2str(k) ')'],...
        ['EEG: F4, LFP: ' lfp_chan],['Baseline: ' num2str(baseline_win)]})
    fnm = ['ETA_' event '_' patnm '_' side '_LFP_' lfp_chan '_WIN' num2str(time([1 end])) 'BAS' num2str(baseline_win) ];
    saveas(fig,fullfile(figdir,[fnm '.jpg']))
    saveas(fig,fullfile(figdir,[fnm '.fig']))
    saveas(fig,fullfile(figdir,[fnm '.pdf']))
    close(fig)
end

end



%--------------------------------------------------------------------------
function eegdats_filt = filtdat(eegdats,freqs,srate)

eegdats_pad = cat(2,eegdats(:,1), eegdats, eegdats(:,end));
eegdats2filt = eegdats_pad(:);
[b,a]=butter(2,[freqs(1)/(srate/2) freqs(2)/(srate/2)],'bandpass');
eegdats_filtpad = filtfilt(b,a,eegdats2filt);
eegdats_fp = reshape(eegdats_filtpad,size(eegdats_pad));
eegdats_filt = eegdats_fp(:,2:end-1);
end