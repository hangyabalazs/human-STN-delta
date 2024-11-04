function EEG_LFP_Wcoh_PD(EventTypes,SubEventTypes)
% EEG_LFP_WCOH_PD  Mean-squared wavelet coherence (MSWC) map between intraop EEG and LFP
%   EEG_LFP_WCOH_PD(EventTypes,SubEventTypes) draws an MSWC map
%   patient-by-patient and averaged across patients around EVENTTYPES. 
%   Compares epochs partitioned based on SUBEVENTTYPES. 
% 
%   Input parameters:
%     EVENTTYPES        1xN cell array of event labels, ex: {'StimulusOn','StopSignal'};
%
%     SUBEVENTTYPES     Nx2 cell array of partition ("subevent") labels, each row
%                       corresponds to an event label, columns to partitions
%                       {'FailedStopTrial','SuccesfulStopTrial';'FailedStopTrial','SuccesfulStopTrial'};
%
% See also: WCOHERENCE, STATCONDFIELDTRIP

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


global filesdir

% Parameters
plot_win = [-1 1];
baseline_win = [-1 -.5];
all_freqs = [1 80];
freqs = [1 4];

timewin_small = [-.5 0];
freq_small = [1.5 2.5];

sess2analyse = getdata2analyse(filesdir, 'rectype','EEG',...
    'rectime','intraop','patients', 'allpatients', 'side','left', 'condition','stimoff');



% WCoherence patient-by-patient - average across event epochs (plot + save struct)
for ei = 1:length(EventTypes)
    
    event = EventTypes{ei}; fprintf('Event: %s,',event);
    
    for sei = 1:3
        if sei<3
            subevent = SubEventTypes{ei,sei}; 
            downsamp = true;
        else
            subevent = '';
            downsamp = false;
        end
        fprintf(' %s...\n',subevent);
        
        wcoh_onebyone(sess2analyse,event,subevent,1,all_freqs,downsamp,0.05,plot_win,baseline_win,1)
        
    end
end



% WCoherence - average across patients (plot + save struct)
for ei = 1:length(EventTypes)
    
    event = EventTypes{ei}; fprintf('Event: %s,',event);

    wcoh_avg(event,'freqs',all_freqs,'alpha',0.05,'plot_win',plot_win,'baseline_win',baseline_win);

end




% Compare partitions


for ei = 1:length(EventTypes)
    event = EventTypes{ei};
    subevs = SubEventTypes(ei,:);
    
    % Compare wavelet coherence map (permutation test with cluster correction)
    compare_wcohmap_patientavg(event,subevs,plot_win,freqs,'all','all',0.05,'cluster',baseline_win)
    
    % Compare selected time windows (paired T test)
    compare_partitions(event,subevs,freq_small,{[]},false,timewin_small,{},baseline_win,'all')
end


end



%--------------------------------------------------------------------------
function wcoh_onebyone(sess2analyse,event,subevent,isfig,freqs,downsamp,alpha,plot_win,baseline_win,chanmean)

global figdir_pd

if isempty(subevent)
    figdir = fullfile(figdir_pd,'introp_EEG_LFP','WCoh',event);
    subevent= event;
else
    figdir = fullfile(figdir_pd,'introp_EEG_LFP','WCoh',[event '_' subevent]);
end
if ~isfolder(figdir); mkdir(figdir); end;

time = -2 :1/250: 2; time = time(1:end-1);
plotlim_inx = dsearchn(time',plot_win'); plotinx = plotlim_inx(1):plotlim_inx(2);
pltime = time(plotinx);
if ~isempty(baseline_win)
    baslim_inx = dsearchn(time',baseline_win'); basinx = baslim_inx(1):baslim_inx(2);
else
    basinx = plotinx;
end

wcoh_structname = fullfile(figdir, ['WavCH_WIN' num2str(plot_win) '_BAS' num2str(baseline_win)]);
stat_structname = fullfile(figdir, ['WavCH_STAT_WIN' num2str(plot_win) '_BAS2' num2str(baseline_win)]);

if exist([wcoh_structname '.mat'])~=2
    WavCH  = struct;
else
    load([wcoh_structname '.mat']);
end

if exist([stat_structname '.mat'])~=2
    WavCH_STAT = struct;
else
    load([stat_structname '.mat']);
end

eeg_chan = 'F4';
for k = 1:length(sess2analyse)
    
    patnm = sess2analyse(k).patient;
    side = sess2analyse(k).side;
    fprintf('%s %s ...\n',patnm,side);
    
    %% EEG
    
    % Load EEG data
    curr_resdir_eeg = sess2analyse(k).folder;
    
    
    [wcoh,phlag,f,coi,chnames,~] = get_patient_wcoh(patnm,side,curr_resdir_eeg,event,subevent,downsamp);
    
    if isempty(wcoh);
        continue;
    end
    
    if chanmean
        wcoh000 = cell(1,length(chnames));
        for ci = 1:length(chnames);
            wcoh000{ci} = cat(3,wcoh{:,ci});
        end
        wcoh_chan = nanmean( cat(4,wcoh000{:}) ,4);
        chnames = {'chanmean'};
    end
    
    ff = f{1}; coicoi = coi{1};
    
    for ci = 1:length(chnames);
        lfp_chan = chnames{ci};
        
        % Average wavelet coherence over epochs
        if ~chanmean
            wcoh_chan = cat(3,wcoh{:,ci});
        end
        
        if ~isempty(baseline_win)
            wcoh_bas_avg = repmat( nanmean(wcoh_chan(:,basinx,:),[2 3]) ,[1 size(wcoh_chan,2) 1]);
            wcoh_bas_std = repmat( std(wcoh_chan(:,basinx,:),[],[2 3],'omitnan') ,[1 size(wcoh_chan,2) 1]);
%                         wcoh_chan = (wcoh_chan - wcoh_bas_avg)./ wcoh_bas_avg;
            
            wcoh_chanavg = (nanmean(wcoh_chan,3)-wcoh_bas_avg)./wcoh_bas_std;
        else
            wcoh_chanavg = nanmean(wcoh_chan,3);
        end
        
%             wcoh_chanavg = nanmean(wcoh_chan,3);
        % Average phase lag
        phlag_chan = cat(3,phlag{:,ci});
        phlag_chanavg = angle(  nanmean( exp(1i*phlag_chan) ,3 )  );
%         ITphlag_chanavg = abs(  nanmean( exp(1i*phlag_chan) ,3 )  );
               
        % Stat
        if ~isnan(alpha)
%              [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(wcoh_chan(:,plotinx,:),ff,alpha,2000,false,'fdr');
            [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(wcoh_chan(:,plotinx,:),ff,alpha,1000,false,'fdr',[],1:length(basinx));
            
            
            % Save  coherence map STATs
            WavCH_STAT.(patnm).(side).([eeg_chan '_' lfp_chan]).mask_ersp = maskersp;
            WavCH_STAT.(patnm).(side).([eeg_chan '_' lfp_chan]).p_ersp = exactp_ersp;
            WavCH_STAT.(patnm).(side).([eeg_chan '_' lfp_chan]).alphafdr = alphafdr;
            WavCH_STAT.(patnm).(side).([eeg_chan '_' lfp_chan]).alpha = alpha;
            
            
            if ~downsamp
                save([stat_structname '.mat'],'WavCH_STAT');
            else
                save([stat_structname '_downsamp.mat'],'WavCH_STAT');
            end
            
            
        end
        
        f_ind = find(ff>freqs(1)&ff<freqs(2));
        fff = ff(f_ind);
        if isfig
            % TIME-FREQ. W-COH PLOT
            wcoh_fig(wcoh_chanavg(f_ind,plotinx),phlag_chanavg(f_ind,plotinx),pltime,fff,coicoi(plotinx));
            fig = gcf;
            if ~isnan(alpha);
                hold on;
                contour(pltime,log2(fff),maskersp(f_ind,:),'Color','red');
            end
            
            % Save plot
            fnm= [patnm '_' event '_' eeg_chan '_' lfp_chan '_F' num2str(freqs) '_WIN' num2str(plot_win) '_BAS2' num2str(baseline_win)];
%             fnm= [patnm '_' event '_' eeg_chan '_' lfp_chan '_F' num2str(freqs) '_WIN' num2str(plot_win) '_BAS' num2str(baseline_win)];
            saveas(gcf,fullfile(figdir,[fnm '.jpg']))
            saveas(gcf,fullfile(figdir,[fnm '.fig']))
            saveas(gcf,fullfile(figdir,[fnm '.pdf']))
            close(gcf);
            
%             % PHASE LAG consistency (resultant vector length)
%             fig = figure;
%             wcoh_fig(ITphlag_chanavg(f_ind,:),phlag_chanavg(f_ind,:),time,fff,coicoi);
%             colormap(pink)
%             caxis([0 0.3])
%             
%             % Save plot
%             fnm= [patnm '_' event '_' eeg_chan '_' lfp_chan '_F' num2str(freqs) '_meanPhLag'];
%             saveas(gcf,fullfile(figdir,[fnm '.jpg']))
%             saveas(gcf,fullfile(figdir,[fnm '.fig']))
%             close(gcf);
        end
        
        
        
        % Save averaged coherence map
        WavCH.(patnm).(side).([eeg_chan '_' lfp_chan]) = wcoh_chanavg;
        
        if ~downsamp
            save([wcoh_structname '.mat'],'WavCH');
        else
            save([wcoh_structname '_downsamp.mat'],'WavCH');
        end
        save(fullfile(figdir,'WavCH_f.mat'),'ff');
        save(fullfile(figdir,'WavCH_coi.mat'),'coicoi');
    end
end
end



%--------------------------------------------------------------------------
function wcoh_fig(wc,phL,t,f,coi)

% fig = figure;
% set(fig,'Visible','on')
h = pcolor(t,log2(f),wc); 
colormap(pink)
h.EdgeColor = 'none';
ax = gca;
ytick=round(pow2(ax.YTick),3);
ax.YTickLabel=ytick;
xlabel('Time(s)')
ylabel('Freq (Hz)')

hcol = colorbar;
hcol.Label.String = 'Magnitude-Squared Coherence';
hold on;
if ~isempty(coi)
    plot(t,log2(coi),'w--','linewidth',2)
end

% if ~isempty(phL)
%     plot_wcoherence_phaselag(phL,wc,3,t,f)
% end

% caxis([0.1 0.4])
caxis([-.4 .4])
end



%--------------------------------------------------------------------------
function [wcoh_pats, wcavg, pats2,fff] = wcoh_avg(event,varargin)

%%
prs = inputParser;
addRequired(prs,'event',@ischar)
addParameter(prs,'subevent','',@ischar)
addParameter(prs,'patients','all',@(x) ischar(x)|iscell(x)) % 'all' | 'tremor-dominant'|'akinetic-rigid'|'mixed' | 'RTdecrease'| 'RTincrease'
addParameter(prs,'region',{},@iscell) % 'all' | {'Motor'} | {'Associative'} | {'Limbic'}
addParameter(prs,'close2centr','all',@ischar) % 'all' | 'Motor' | 'Associative' | 'Limbic'
addParameter(prs,'maxpower','all',@ischar) % 'all' | 'beta' | 'high_delta' 
addParameter(prs,'chanmean',1,@isnumeric)
addParameter(prs,'side','left',@(x) ischar(x)|iscell(x)) % 'both' | {'left','right'}
addParameter(prs,'iscell',true,@(x) islogical(x)| isnumeric(x)) % [] | true | false
addParameter(prs,'freqs',[1 100],@(x) isvector(x)|isnumeric(x));
addParameter(prs,'isfig',true,@islogical);
addParameter(prs,'downsamp',false,@islogical);
addParameter(prs,'alpha',NaN,@isnumeric);
addParameter(prs,'plot_win',[-2 2],@isvector);
addParameter(prs,'baseline_win',[],@(x) isvector(x)|isnumeric(x));
parse(prs,event,varargin{:})
g = prs.Results;

global rootdir figdir_pd

eegchan = 'F4';
time = -2 :1/250: 2; time = time(1:end-1);
pltime = g.plot_win(1):1/250:g.plot_win(2);
timinx = dsearchn(time',pltime');
if ~isempty(g.baseline_win)
    bastime = g.baseline_win(1):1/250:g.baseline_win(2);
    basinx = dsearchn(time',bastime');
    basinx = basinx- timinx(1)+1;
end

if isempty(g.subevent)
    figdir = fullfile(figdir_pd,'introp_EEG_LFP','WCoh',event);
else
    figdir = fullfile(figdir_pd,'introp_EEG_LFP','WCoh',[event '_' g.subevent]);
end

if ~isfolder(figdir); mkdir(figdir); end;

wcoh_structname = fullfile(figdir, ['WavCH_WIN' num2str(g.plot_win) '_BAS']);

% wcoh_structname = fullfile(figdir, ['WavCH_WIN' num2str(g.plot_win) '_BAS' num2str(g.baseline_win)]);

if g.downsamp
    load([wcoh_structname '_downsamp.mat']);
else
    load([wcoh_structname '.mat']);
end

load(fullfile(figdir,'WavCH_f.mat'));
load(fullfile(figdir,'WavCH_coi.mat'));


% wcoh_avg_structname = fullfile(figdir, ['WavCH_avg_WIN' num2str(g.plot_win) '_BAS' num2str(g.baseline_win) '.mat']);
% % wcoh_avg_structname = fullfile(figdir, ['WavCH_avg_WIN' num2str(g.plot_win) '_BAS.mat']);
% 
% if exist(wcoh_avg_structname)~=2
%     WavCH_avg  = struct;
% else
%     load(wcoh_avg_structname);
% end


% Freqs
if ~isempty(g.freqs)
    f_inx = ff>g.freqs(1)&ff<g.freqs(2);
    fff = ff(f_inx);
else
    f_inx = 1:length(ff);
    fff = ff;
end

if strcmp(g.patients,'all')
    pats = fieldnames(WavCH);
else
    
    goodpats = fieldnames(WavCH);
    grpats = clinical_groups({g.patients});
    grpats = cat(2,grpats{:});
    pats = intersect(goodpats,grpats); % some patients might have missing data
end

if strcmp(g.side,'both')
    sids = {'left','right'};
else
    if ~iscell(g.side); sids = {g.side}; else; sids = g.side; end;
end

k = 1; wcoh = {}; pats2 = {};
for pp = 1:length(pats)
    patnm = pats{pp};
    for ss = 1:length(sids)
        sidenm = sids{ss};
        
        %% Select channels to average
        inchans = {};
        subreg_tit = '';
%         % Get channels in predefined subregion
%         if ~isempty(g.region{1});
%             inchans0 = cell(1,length(g.region));
%             for rr = 1:length(g.region)
%                 try
%                     inchans0{rr} = get_chan_subregion(patnm,sidenm,g.region{rr});
%                 catch
%                     disp('')
%                 end
%             end
%             inchans = cat(2,inchans0{:});
%             inchans1 = inchans;
%             subreg_tit = [g.region{:}];
%         end
        
        
        if ~g.chanmean
            % Get channel closest to centroid of one STN subreg.
            if  ~strcmp(g.close2centr,'all')
                load(fullfile(rootdir,'Channels_close2centroids.mat'));
                inchans = centr_tab{[patnm '_' sidenm(1)],g.close2centr};
                subreg_tit = ['Close to ' g.close2centr ' centr.'];
                
            end
            
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
                    [~,~,cell_chan(cc),~] = cellid2tags(cellids{cc});
                end
                inchans = arrayfun(@(x) ['Ch' num2str(x)],unique(cell_chan),'UniformOutput',0);
                
                if ~g.iscell
                    fdnm = fieldnames(WavCH.(patnm).(sidenm));
                    allchans = cellfun(@(x) x(strfind(x,'_')+1:end),fdnm,'UniformOutput',0);
                    
                    inchans00 = inchans;
                    inchans = allchans( ~ismember(allchans,inchans00) );
                end
                inchans2 = inchans;
                cell_tit = ['cells-' char(string(g.iscell))];
            else
                cell_tit = '';
            end
            
            % If no condition, get all channels
            if ~isempty(g.region{1})&& ~isempty(g.iscell)
                inchans = intersect(inchans1,inchans2);
                
                %         elseif isempty(g.region{1})&& isempty(g.iscell) && strcmp(g.close2centr,'all') && strcmp(g.maxpower,'all');
                %             fdnm = fieldnames(WavCH.(patnm).(sidenm));
                %             inchans = cellfun(@(x) x(strfind(x,'_')+1:end),fdnm,'UniformOutput',0);
            end
        else
            %             inchans = {'chanmean'};
             fdnm = fieldnames(WavCH.(patnm).(sidenm));
             inchans = cellfun(@(x) x(strfind(x,'_')+1:end),fdnm,'UniformOutput',0);

        end
        
        
        
        wcoh_chans = cell(1,length(inchans));
        for ci = 1:length(inchans)
            lfpchan = inchans{ci};
            if isfield(WavCH.(patnm).(sidenm),[eegchan '_' lfpchan])
                wcch = double(WavCH.(patnm).(sidenm).([eegchan '_' lfpchan]));
                
                if isempty(g.baseline_win)
                    % Remove outliers
                    wcch(abs(wcch)>1) = NaN;
                    wcch(isnan(wcch)) = nanmean(wcch(:));
                end
                
                wcoh_chans{ci} = wcch;
            else
                fprintf('No %s in %s %s\n',lfpchan,patnm,sidenm);
            end
        end
        if ~isempty(cat(3,wcoh_chans{:}))
            wcoh{k} =  mean(cat(3,wcoh_chans{:}),3);
            k = k+1;
            pats2 = cat(2,pats2,patnm);
        end
    end
end
wcoh_pats = cat(3,wcoh{:});

wcoh_pats = wcoh_pats(f_inx,timinx,:);

if ~isnan(g.alpha)
    
    [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(wcoh_pats,fff,g.alpha,1000,false,'fdr',[],1:length(basinx));
    
end

if ~isempty(g.baseline_win)
    bas_avg = repmat(mean(wcoh_pats(:,basinx,:),[2 3]), [1 size(wcoh_pats,2) 1]);
    bas_std = repmat(std(wcoh_pats(:,basinx,:),[], [2 3],'omitnan'), [1 size(wcoh_pats,2) 1]);
    wcavg = (mean(wcoh_pats,3) - bas_avg)./bas_std;
else
    wcavg = mean(wcoh_pats,3);
end



nnr = size(wcoh_pats,3);
% WavCH_avg.(event).(['F' num2str(g.freqs)]).([subreg_tit '_' cell_tit]) = wcoh_avg;
% save(fullfile(figdir,'WavCH_avg.mat'),'WavCH_avg');
%%
if g.isfig
    wcoh_fig(wcavg,[],pltime,fff,coicoi(timinx));
    fig = gcf;
    if ~isnan(g.alpha)
        hold on;
        contour(pltime, log2(fff),maskersp,'Color','red');
    end
    
    if ~g.chanmean
        
        title({[g.patients ' patients (n=' num2str(nnr) ')'],event,['EEG: F4, LFP: ' subreg_tit ', ' cell_tit],[num2str(g.freqs) ' Hz']})
        fnm= ['AVG_' event '_' subreg_tit '_' cell_tit '_F' num2str(g.freqs) '_PATS_' g.patients '_DS' char(string(g.downsamp)) '_WIN' num2str(g.plot_win) '_BAS' num2str(g.baseline_win)];
    else
        title({[g.patients ' patients (n=' num2str(nnr) ')'],event,['EEG: F4, LFP: chanmean'],[num2str(g.freqs) ' Hz']})
        fnm= ['AVG_' event '_chanmean_F' num2str(g.freqs) '_PATS_' g.patients '_DS' char(string(g.downsamp)) '_WIN' num2str(g.plot_win) '_BAS' num2str(g.baseline_win)];
    end
    saveas(gcf,fullfile(figdir,[fnm '.jpg']))
    saveas(gcf,fullfile(figdir,[fnm '.fig']))
    saveas(gcf,fullfile(figdir,[fnm '.emf']))
    close(gcf);
end
end



%--------------------------------------------------------------------------
function compare_partitions(event,subevs,freqs,region,iscell,timewin,patgroup_nm,baseline_win,maxpower)


global figdir_pd

twin_nr = size(timewin,1);
part_nr = 2;
srate = 250;

% Indices of data points corresponding to predef timewindow
plot_win = [-1 1];
alltime = plot_win(1):1/srate:plot_win(2);
[timinx,timelab]= deal(cell(1,twin_nr));
for t = 1:twin_nr
    tinx = dsearchn(alltime',timewin(t,:)');
    timinx{t} = tinx(1):tinx(2);
    timelab{t} = [num2str(timewin(t,1)) ' - ' num2str(timewin(t,2))];
end


%%
[wcoh_patsall,allpats,wcoh_pats,pats] = deal(cell(1,part_nr));
for rr = 1:part_nr
    
    [wcoh_patsall{rr}, ~,allpats{rr}] = wcoh_avg(event,'subevent',subevs{rr},...
        'region',region,'iscell',iscell,'freqs',freqs,'isfig',false,'downsamp',true,...
        'baseline_win',baseline_win,'plot_win',plot_win,'maxpower',maxpower);
end

% Discard patients with trials in only 1 contingency
[pats0,pix1,pix2] = intersect(allpats{:});
wcoh_pats{1} = wcoh_patsall{1}(:,:,pix1);
wcoh_pats{2} = wcoh_patsall{2}(:,:,pix2);
pats{1} = pats0; pats{2} = pats0;
%%
% Create cell array for grouped boxplots + ANOVA
bp = cell(1,twin_nr);
for tk = 1:twin_nr
    % Average over time (one timewindow) and frequency components
    bp{tk} = cellfun(@(x) squeeze(mean(x(:,timinx{tk},:),[1 2])),wcoh_pats,'UniformOutput',0);
    
    % Complete matrices with NaN if there are missing values
    mL = max(cellfun(@length,bp{tk}));
    bp2{tk} = cell2mat(cellfun(@(x) cat(1,x, nan(mL-length(x),1) ),bp{tk},'UniformOutput',0));
    
    % TTEST (compare partitions from one time-window)
    [~,pval_parts(tk),~,~] = ttest2(bp2{tk}(:,1),bp2{tk}(:,2));
end

bp3 = bp2{1};
fig = figure;
% boxplotGroup(bp2,'secondaryLabels',subevs,'primaryLabels', timelab);
boxplot(bp3,subevs); hold on;
for kk = 1:size(bp3,2);
    for kkk = 1:size(bp3,1)
        scatter(kk+randi([-100 100],1)*0.0015,bp3(kkk,kk),[],'k','filled'); hold on;
    end
end
set_my_boxplot(gca)
% yL = [0 .8]; %[0.2 0.4];
% ylim(yL)
yL = ylim;

if twin_nr>1
    % TTEST (compare time windows for each subregion)
    for jk = 1:part_nr
        
%         [~,pval_ttest(jk),~,~] = ttest2(bp2{1}(:,jk),bp2{2}(:,jk));
        [~,pval_ttest(jk),~,~] = ttest(bp2{1}(:,jk),bp2{2}(:,jk));
        
        % Write ttest p to plot
        if pval_ttest(jk)<0.05; col = 'r'; else; col = 'k'; end;
        text((jk-1)*2+jk,yL(2)*0.9,['T1 vs T2: p=' num2str(pval_ttest(jk))],'Color',col)
    end
else
    % Lines for patients
    line(repmat([1; 3],[1 size(bp3,1)]),bp3','Color',[.5 .5 .5]);
end

% Write Anova p to plot
for tk = 1:twin_nr
    
    if pval_parts(tk)<0.05; col = 'r'; else; col = 'k'; end;
    text(1,yL(1)+yL(1)*tk*0.1,['T' num2str(tk) ': ' subevs{1}(1:4) ' vs ' subevs{2}(1:4) ' p=' num2str(pval_parts(tk))],'Color',col)
end


% Mark patient groups
patgroup_lab = '';
if ~isempty(patgroup_nm)
    hold on;
    patgroups = clinical_groups(patgroup_nm);
    cols = jet(length(patgroups)+2);
    Len = size(bp2{1},1);
    for rr = 1:part_nr
        pats_grouped = cellfun(@(x) find(ismember(pats{rr},x)),patgroups,'UniformOutput',false);
        
        bp3 = cell2mat(cellfun(@(x) x(:,rr) ,bp2,'UniformOutput',0));
        %         bp4 = cat(2,nan(Len,(rr-1)*3),bp3);
        for gr = 1:length(pats_grouped)
            % dat = bp4(pats_grouped{gr},:);
            % sc{gr} = datpointsplot(dat,(rr-1)*3+[1,2],{},cols(gr,:));
            
            dat = bp3(pats_grouped{gr},:);
            if twin_nr>1
                sc{gr} = line(repmat((rr-1)*3+[1,2],[size(dat,1),1])',dat','Color',cols(gr,:),'LineWidth',2)
            else
                sc{gr} = scatter(ones(1,size(dat,1))*(rr-1)*2+1+randi([-1 1],1,size(dat,1))*0.1,dat,[],cols(gr,:),'filled'); hold on;
                
            end
            
        end
        ise = cellfun(@isempty,sc);
        scleg = cellfun(@(x) x(1),sc(~ise));
        if rr==1
            legend(scleg,patgroup_nm(~ise),'AutoUpdate','off')
        end
        
        
    end
    
    patgroup_lab = cellfun(@(x) x(1:2),patgroup_nm,'UniformOutput',0);
    patgroup_lab = [patgroup_lab{:}];
end

if ~isempty(iscell)
    cell_tit = ['cells-' char(string(iscell))];
else;
    cell_tit = '';
end
title(timelab)

figdir_stat = fullfile(figdir_pd,'introp_EEG_LFP','WCoh',event,'Stat');
if ~isfolder(figdir_stat); mkdir(figdir_stat); end;

partilab = cellfun(@(x) x(1),subevs,'UniformOutput',0); partilab = [partilab{:}];

fnm = [partilab '_PART_'  event '_TIME' timelab{:} '_FREQ' num2str(freqs) '_' cell_tit '_' patgroup_lab '_BAS' num2str(baseline_win) '_MAXPOW' maxpower];
saveas(fig,fullfile(figdir_stat,[fnm '.jpg']))
saveas(fig,fullfile(figdir_stat,[fnm '.fig']))
saveas(fig,fullfile(figdir_stat,[fnm '.pdf']))
close(fig)



end



%--------------------------------------------------------------------------
function compare_wcohmap_patientavg(event,subevs,plot_win,freqs,close2centr,maxpower,alpha,mcorrect,baseline_win)

global figdir_pd
figdir = fullfile(figdir_pd,'introp_EEG_LFP','WCoh',event);
figdir_stat = fullfile(figdir,'Stat');

load(fullfile(figdir,'WavCH_coi.mat'));


% Time
sr = 250;

pltime = plot_win(1):1/sr:plot_win(2);
basinx = dsearchn(pltime',baseline_win');
%%

[basdat, ~,allpats_b,f] = wcoh_avg(event,'subevent','',...
    'region',{[]},'iscell',[],'freqs',freqs,'isfig',false,'downsamp',false,...
    'close2centr',close2centr,'maxpower',maxpower,'baseline_win',[],'plot_win',plot_win);
basavg =  repmat(nanmean(basdat(:,basinx,:),2),[1 size(basdat,2) 1]);
bassd =  repmat(std(basdat(:,basinx,:),[],2,'omitnan'),[1 size(basdat,2) 1]);


statdat = {};
for ss = 1:2
    
    [sdat, ~,allpats{ss},ff] = wcoh_avg(event,'subevent',subevs{ss},...
        'region',{[]},'iscell',[],'freqs',freqs,'isfig',false,'downsamp',true,...
        'close2centr',close2centr,'maxpower',maxpower,'baseline_win',[],'plot_win',plot_win);
    
    patix = ismember(allpats_b,allpats{ss});
    
    sdat2 = (sdat-basavg(:,:,patix))./ bassd(:,:,patix);
    statdat{ss,1} = sdat2;
    
    sdat = [];
end




                
 
% 
% [pcond, ~, ~, statscond, ~, ~] = std_stat(statdat,'condstats','on','groupstats','off','mode','eeglab',...
%     'method','perm','mcorrect',mcorrect,'alpha',alpha);



[F, df, pval]  = statcondfieldtrip(statdat, 'method','montecarlo','naccu', 1000,...
    'mcorrect',mcorrect,'alpha',alpha);
pcond{1} = pval<alpha;
statscond{1} = F;
                

fig = figure;
for sei = 1:2
    subplot(3,1,sei)
    wcoh_fig(mean(statdat{sei,1},3),[],pltime,ff,[]);
caxis([-20 20])
    title(subevs{sei})
end
subplot(3,1,3)
wcoh_fig(statscond{1},[],pltime,ff,[]); caxis([-2 2])
title('F-S: T values')

hold on; contour(pltime,log2(ff),pcond{1},'Color','red');

fnm = ['mapFS2_PART_' event '_TIME' num2str(plot_win) '_FREQ' num2str(freqs) '_CLOSE2_' close2centr '_MAXPOW' maxpower '_ALPHA' num2str(alpha) '_' mcorrect '_BAS' num2str(baseline_win)];
saveas(fig,fullfile(figdir_stat,[fnm '.jpg']))
saveas(fig,fullfile(figdir_stat,[fnm '.fig']))
saveas(fig,fullfile(figdir_stat,[fnm '.pdf']))
close(fig)

% 
% fig = figure;
% wcoh_fig(statscond{1},[],pltime,ff,coicoi(timinx)); caxis([-2 2])
% colormap(pink)
% title('F-S: T values')
% colormap(pink)
% hold on; contour(pltime,log2(ff),pcond{1},'Color','red');
% 
% 
% fnm = ['mapFS_PART_' event '_TIME' num2str(plot_win) '_FREQ' num2str(freqs) '_DIFFMAP'];
% saveas(fig,fullfile(figdir_stat,[fnm '.jpg']))
% saveas(fig,fullfile(figdir_stat,[fnm '.fig']))
% saveas(fig,fullfile(figdir_stat,[fnm '.eps']))
% close(fig)

end



