function spike_triggered_average_PD(sess2analyse,chanmean,freq,save_window,plot_window,EventTypes)
%SPIKE_TRIGGERED_AVERAGE_PD  Spike-triggered average
%  SPIKE_TRIGGERED_AVERAGE_PD(sess2analyse,chanmean,freq,save_window,plot_window,EventTypes)
%       -Saves spike-triggered LFP epochs for sessions specified in SESS2ANALYSE, 
%       in SAVE_WINDOW time window
%       -Plots spike-triggered LFP average with standard error in PLOT_WINDOW time window, including
%       units significantly coupled to FREQ Hz (calcaulted around
%       EVETTYPES).
%
% Required inputs:
%     SESS2ANALYSE  struct containing all necessary information (name of patient, side
%                    of experiment, tag of condition, session folder path) of
%                    sessions that need to be analysed (see GETDATA2ANALYSE)
%     CHANMEAN      if true, channel averaged LFP data is used, if false
%                   LFP data is plotted channel-by-channel
%                   (relevant only for intraop LFP data)
%                   true | false 
%     FREQ          1x2 vector, boundaries of frequency band of
%                   interest (Hz) (first column is
%                   the lower limit, second column is the lower limit)
%     SAVE_WINDOW   1x2 vector, time window limit (sec) around spikes to
%                   save spike-triggered epochs
%     PLOT_WINDOW   1x2 vector, time window limit (sec) around spikes to
%                   plot spike-triggered average
%     EVENTTYPES    cell array of event labels, to select significantly
%                   coupled units around specified events
%
% See also: ASTRANORM2, PC_CELL_LEVEL, ERRORSHADE

% Johanna Petra Szabó, Hangya Balázs, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


if chanmean; chm = 1; else; chm =0; end;

% Save spike-triggered epochs
STA_PD(sess2analyse,false,chm,freq,save_window)


% Plot spike-triggered averages for significantly coupled units
rectype = sess2analyse(1).rectype;
STavg_signPC(rectype,chm,EventTypes,freq,plot_window)


end





function STA_PD(sess2analyse,isplot,chanmean,freq,save_window);


global cell_dir

rectype = sess2analyse(1).rectype;
figdir = fullfile(cell_dir,[rectype '_STA_STS'],'STA');
fastif(~isdir(figdir),mkdir(figdir),0);

for snr = 1:length(sess2analyse)
    
    curr_resdir = sess2analyse(snr).folder;
    side = sess2analyse(snr).side;
    condition = sess2analyse(snr).tag;
    patnm = sess2analyse(snr).patient;
    currsess = sess2analyse(snr).sessfolder;
    
    
    sL = strfind(currsess,filesep); sessnm = currsess(sL(end)+1:end);
    cellids = findcell('rat',patnm,'session',sessnm);
    
    
    fprintf('%s %s %s...\n',patnm, side, condition);
    
    % load LFP/EEG
    eegfile = dir(fullfile(curr_resdir,'EEG_*_05HP.set'));
    if isempty(eegfile)
        fprintf('No EEG file %s %s\n',patnm, sessnm)
        continue
    end
    EEG = pop_loadset(fullfile(eegfile.folder,eegfile.name));
    
    if ~isempty(freq)
        [EEG, com, b] = pop_eegfiltnew(EEG,freq(1),freq(2));
    end
    
    eeg_srate = EEG.srate;
    eeg_times = EEG.xmin:1/eeg_srate:EEG.xmax;
    
    
    cells = findcell('rat',patnm,'session',sessnm);
    
    
    for iC = 1:length(cells)
        
        
        cellid = cells{iC};
        fprintf('CellID: %s...\n',cellid)
        
        % get LFP/EEG subdirectory
        [~,~,channr,unit] = cellid2tags(cellid);
        
        % load spike timestamps
        
        load(fullfile(currsess,['Ch' num2str(channr) '_' num2str(unit)]));
        
        % load matching channel data (if not channel mean is used)
        if strcmp(rectype,'LFP')
            chlab = ['chm' num2str(chanmean)];
            if chanmean~=1
                try
                    chinx = find(strcmp(['Ch' num2str(channr)],{EEG.chanlocs(:).labels}));
                catch
                    chinx = channr;
                end
                dat = EEG.data(chinx,:);
            else
                dat = mean(EEG.data,1);
            end
        elseif strcmp(rectype,'EEG')
            chlab = EEG.chanlocs.labels;
            dat = EEG.data(1,:,:);
        end
        
        spike_inx = dsearchn(eeg_times',TS/10000);
        
        wn = diff(save_window);
        
        
        [sta, ~, ~, ~, ~, ~, st] = astanorm2(spike_inx',dat,eeg_srate,wn*eeg_srate);
        
        win_times = save_window(1):1/eeg_srate:save_window(2);
        
        % Standard error
        sta_sd = std(st, [], 1);
        sta_sd = sta_sd / sqrt(size(st, 1));
        
        if isplot
            % eta figure
            fig = figure;
            errorshade(win_times,sta,sta_sd,'LineColor',[0 0 0.6],'ShadeColor',[0 0 0.6],'LineWidth',1.5)
            
            
            
            cellidtit = cellid; cellidtit(ismember(cellid,'_')) = '-';
            title(cellidtit)
            set(fig,'Position',get(0,'Screensize'));
            
            saveas(fig,fullfile(figdir,[cellid '_STA_' chlab '_FR' num2str(freq) '.jpg']));
            savefig(fig,fullfile(figdir,[cellid '_STA_' chlab '_FR' num2str(freq) '.fig']));
            close(fig)
        end
        
        if exist(fullfile(figdir,['STAmat_' chlab '_FR' num2str(freq) '.mat']),'file')~=2
            STAmat = struct;
        else
            load(fullfile(figdir,['STAmat_' chlab '_FR' num2str(freq) '.mat']));
        end
        
        cellidF = cellid; cellidF(strfind(cellid,'.')) = '_';
        STAmat.(cellidF).times = win_times;
        STAmat.(cellidF).sta = sta;
        STAmat.(cellidF).sta_sd = sta_sd;
        save(fullfile(figdir,['STAmat_' chlab '_FR' num2str(freq) '.mat']),'STAmat')
        
    end
end
end


%--------------------------------------------------------------------------
function STavg_signPC(rectype,chanmean,EventTypes,freq,plot_window)

global cell_dir
PCdir = fullfile(cell_dir,'PC',[rectype '_phase_coupling']);
frnm = 'dom_high_delta';
PC_win = [-1.5 1.5];

if contains(rectype,'LFP')&& chanmean==1
    chtit = ['chanmean_allSTN'];
elseif contains(rectype,'LFP')&& chanmean==0
    chtit = 'by-channel';
elseif strcmp(rectype,'EEG')
    chtit = 'F4';
end

for ei = 1:length(EventTypes)
    event = EventTypes{ei};
    
    resdir = fullfile(PCdir,chtit,event,event);
    load(fullfile(resdir,[frnm '_' num2str(PC_win)],'PC_results_dsno_1win.mat'));
    
    cellids0 = fieldnames(PC_results.Hilb_PC);
    ray_sign = structfun(@(x) x.Ray_P<=0.05,PC_results.Hilb_PC);
    cellids = cellids0(ray_sign);
    
    chlab =['chm' num2str(chanmean)];
    fnm = ['signPC_' event '_' chlab '_FR' num2str(freq) '_WIN' num2str(plot_window) '_SE'];
    STA_avg(cellids,rectype,fnm,chlab,freq, plot_window);
end
end

%--------------------------------------------------------------------------

function STA_avg(cellids,rectype,fnm,chlab,freq,plot_window,errshade)

global cell_dir

figdir = fullfile(cell_dir,[rectype '_STA_STS'],'STA');

load(fullfile(figdir,['STAmat_' chlab '_FR' num2str(freq) '.mat']))

eeg_srate = 250;

[sta,sta_sd] = deal([cell(length(cellids),1)]);
for ic = 1:length(cellids)
    cellid = cellids{ic};
    
    cellidF = cellid; cellidF(strfind(cellid,'.')) = '_';
    
    if ~isfield(STAmat,cellidF)
        fprintf('No STA %s\n',cellidF); continue;
    else
        times = STAmat.(cellidF).times;
        newtimes = plot_window(1):1/eeg_srate:plot_window(2);
        newtiminx = dsearchn(times',newtimes');
        
        sta{ic} = double(STAmat.(cellidF).sta(newtiminx));
        sta_sd{ic} = double(STAmat.(cellidF).sta_sd(newtiminx));
    end
end
fig = figure;
stam = cat(1,sta{:});
sta_avg = mean(stam,1);
sta_sd = std(stam,[],1);
sta_se = sta_sd./sqrt(size(stam,1));
errorshade(times(newtiminx),sta_avg,sta_se,'LineColor',[0 0 0.6],'ShadeColor',[0 0 0.6],'LineWidth',1.5)
legend({'mean','SE'})
setmyplot_balazs(gca)
resdir = fullfile(figdir,'STAavg'); if ~isfolder(resdir); mkdir(resdir); end;
saveas(fig,fullfile(resdir,[fnm '.jpg']))
saveas(fig,fullfile(resdir,[fnm '.fig']))
saveas(fig,fullfile(resdir,[fnm '.pdf']))
close(fig);
end