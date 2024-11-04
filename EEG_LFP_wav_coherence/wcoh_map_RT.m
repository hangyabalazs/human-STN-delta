function wcoh_map_RT
%WCOH_MAP_RT   Correlation map between reaction time and wavelet coherence map
%   WCOH_MAP_RT calculates trial-by-trial correlation (Pearson's correlation) 
%   between reaction time and wavelet coherence map of each patient.
%   Significant correlation tested by TTEST on the level of individual patients
%   is marked by contouring on individual patient figures. Significant
%   change in correlation tested by permutation test wih fdr correction is
%   marked by contouring on patient-averaged figures.


% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

global filesdir

plot_window = [-1 1];
baseline_win = [-1 -.5];
event = 'StimulusOn';
freq_bands = [1 80];


sess2analyse = getdata2analyse(filesdir, 'rectype','EEG',...
    'rectime','intraop','patients', 'allpatients', 'side','left', 'condition','stimoff');


% Calculate and save correlation maps patient-by-patient
WCohmap_RTcorr_pats(sess2analyse,event,baseline_win,plot_window)

% Draw correlation maps patient-by-patient
RTcorrmapW_FIG_pats(sess2analyse,event,event,plot_window,freq_bands,'basnorm')

% Draw correlation maps averaged across patients + perform statistics
WCohmap_patavg(sess2analyse,event,plot_window,baseline_win,freq_bands,0.05)
end



%--------------------------------------------------------------------------
function WCohmap_patavg(sess2analyse,event,plot_window,baseline_win,freq_bands,alpha)

global rootdir figdir_pd

figdir = fullfile(figdir_pd,'introp_EEG_LFP','behav_corr','RT_corrmap');
load(fullfile(figdir,'RTcorrMAP.mat'))

epoch_win = [-2 2];
sr = 250;
times = epoch_win(1):1/sr:epoch_win(2);times = times(1:end-1);
newtimes = plot_window(1):1/sr:plot_window(2);newtimes = newtimes(1:end-1);
timinx = dsearchn(times',newtimes');

baslim_inx = dsearchn(newtimes',baseline_win');
basinx = baslim_inx(1):baslim_inx(2);


wdir = fullfile(figdir_pd,'introp_EEG_LFP','WCoh',event);
load(fullfile(wdir,'WavCH_f.mat'));

f_ind = intersect(find(ff>=freq_bands(1)),find(ff<=freq_bands(2)));
fff = ff(f_ind);


sessnr = length(sess2analyse);
for si = 1:sessnr
    
    condi = sess2analyse(si).tag;
    patnm = sess2analyse(si).patient;
    side = sess2analyse(si).side;
    
    coefD{si} = RTcorrMAP.([patnm '_' side]).chanmean.(event).(event).trialnorm.CoefMap(f_ind,timinx,:);
    
    
end
allcoef_stat = cat(3,coefD{:});
%
formula = 'mean(arg1,3);';
[exactp_ersp,maskersp,~] = boostat_eeglab_J(allcoef_stat,[],alpha,1000,...
    false,'fdr',formula,basinx);



fig = figure;
imagesc(newtimes,1:length(fff),mean(allcoef_stat,3));
% crange = caxis;
crange = [-.1 .1];
hold on;
contour(newtimes,1:length(fff),maskersp,'Color','r')
caxis(crange)
yti = round(linspace(1,length(fff), 4)) ;
yticks(yti); yticklabels( arrayfun(@num2str,fff(yti),'UniformOutput',0) );

colormap(bone);
title({['Patietn AVG, channelmean'],[side ', ' condi]});
xlabel('Time (s)'); ylabel('Frequency (Hz)');
c=colorbar; ylabel(c,'R','Rotation',270);
fnm = fullfile(figdir,['AVG_' side '_' condi '_' event '_FR' num2str(freq_bands)]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);

end


%--------------------------------------------------------------------------
function WCohmap_RTcorr_pats(sess2analyse,event,baseline_win,plot_window)

global figdir_pd


figdir = fullfile(figdir_pd,'introp_EEG_LFP','behav_corr','RT_corrmap');
if ~isfolder(figdir); mkdir(figdir); end;

epoch_win = [-2 2];
sr = 250;
times = epoch_win(1):1/sr:epoch_win(2);times = times(1:end-1);
newtimes = plot_window(1):1/sr:plot_window(2); newtiminx = dsearchn(times',newtimes');
baslim_inx = dsearchn(times',baseline_win');
basinx = baslim_inx(1):baslim_inx(2);


sessnr = length(sess2analyse);
for si = 1:sessnr
    
    if exist(fullfile(figdir,'RTcorrMAP.mat'))==2
        load(fullfile(figdir,'RTcorrMAP.mat'))
    else
        RTcorrMAP = struct;
    end
    
    curr_resdir_eeg = sess2analyse(si).folder;
    condi = sess2analyse(si).tag;
    patnm = sess2analyse(si).patient;
    side = sess2analyse(si).side;
    
    fprintf('%s, %s, %s...',patnm, condi, side)
    
    
    
   
    
    
    % MSWC
    [wcoh_alltrials,~,f_c,coi,chnames,TEinx_eeg2] = get_patient_wcoh(patnm,side,...
        curr_resdir_eeg,event,event,false);
    
    wdat0 =  {};
    for ccc = 1:size(wcoh_alltrials,2)
        wdat0{ccc} = cat(3,wcoh_alltrials{:,ccc});
    end
    wdat = nanmean( cat(4,wdat0{:}) ,4);
    
    
     % Reaction Time
    
    [RT,Go_RT,FAlarm_RT,Hit,NoStopTrials StopTrials,perf,stopperf] = ...
        calc_RT(sess2analyse(si),['TrialEvents_' condi '.mat'],false);
    RT2 = RT{1}(TEinx_eeg2);
    
    % MSWC NORM + Correlation
    % Full-length trial norm
    wbasT = nanmean(wdat ,2);
    wsdT = nanstd( wdat ,[],2);
    
    repbas = repmat(wbasT,[1 size(wdat,2)]);
    repsd = repmat(wsdT,[1 size(wdat,2)]);
    
    wdat_Tnorm = (wdat - repbas)./repsd;
    
    RTcorrMAP = RTmap_save(RT2,wdat_Tnorm,'trialnorm',patnm,side,event,event,'chanmean',figdir,RTcorrMAP);

    
    % Baseline norm
%     wbas = mean( nanmean( wdat(:,basinx,:) ,3) ,2);
%     wsd = std(  nanmean(wdat(:,basinx,:) ,3) ,[],2);
%     
%     repbas = repmat(wbas,[1 size(wdat,2:3)]);
%     repsd= repmat(wsd,[1 size(wdat,2:3)]);
    
    wbas = nanmean( wdat(:,basinx,:) ,2);
    wsd = nanstd(  wdat(:,basinx,:) ,[],2);
    
    repbas = repmat(wbas,[1 size(wdat,2)]);
    repsd = repmat(wsd,[1 size(wdat,2)]);
    
    wdat_norm = (wdat - repbas)./repsd;
    
    
    RTcorrMAP = RTmap_save(RT2,wdat_norm,'basnorm',patnm,side,event,event,'chanmean',figdir,RTcorrMAP);
    
end
end



%--------------------------------------------------------------------------
function RTcorrMAP = RTmap_save(RT2,tf_1chan,normtype,patnm,side,event,evty,act_chan,figdir,RTcorrMAP)
rtnn = ~isnan(RT2);

tf = permute(tf_1chan(:,:,rtnn),[3 1 2]);
rtmat = repmat(RT2(rtnn)', [1 size(tf_1chan,[1 2])] );

coefD = []; pvalD =[];
for j = 1:size(tf_1chan,2)
    [coef, pval] = corr(tf(:,:,j),rtmat(:,:,j));
    
    coefD = [coefD diag(coef)];
    pvalD = [pvalD diag(pval)];
end

RTcorrMAP.([patnm '_' side]).(act_chan).(event).(evty).(normtype).CoefMap = coefD;
RTcorrMAP.([patnm '_' side]).(act_chan).(event).(evty).(normtype).PMap = pvalD;

save(fullfile(figdir,'RTcorrMAP.mat'),'RTcorrMAP');

end

%--------------------------------------------------------------------------
function RTcorrmapW_FIG_pats(sess2analyse,event,evty,plot_window,freq_bands,normtype)

global figdir_pd


figdir = fullfile(figdir_pd,'introp_EEG_LFP','behav_corr','RT_corrmap');

wdir = fullfile(figdir_pd,'introp_EEG_LFP','WCoh',event);
load(fullfile(wdir,'WavCH_f.mat'));

f_ind = intersect(find(ff>=freq_bands(1)),find(ff<=freq_bands(2)));



epoch_win = [-2 2];
new_sr = 250;
times = epoch_win(1):1/new_sr:epoch_win(2);times = times(1:end-1);
newtimes = plot_window(1):1/new_sr:plot_window(2);newtimes = newtimes(1:end-1);
timinx = dsearchn(times',newtimes');


load(fullfile(figdir,'RTcorrMAP.mat'))

sessnr = length(sess2analyse);
for si = 1:sessnr
    
    
    condi = sess2analyse(si).tag;
    patnm = sess2analyse(si).patient;
    side = sess2analyse(si).side;
    
    fprintf('%s, %s, %s...',patnm, condi, side)
    
    
    coefD = RTcorrMAP.([patnm '_' side]).chanmean.(event).(evty).(normtype).CoefMap(f_ind,:);
    pvalD = RTcorrMAP.([patnm '_' side]).chanmean.(event).(evty).(normtype).PMap(f_ind,:);
    
    corrmapfig(coefD,pvalD,normtype,times,timinx,ff(f_ind),patnm,condi,side,event,evty,'chanmean',figdir);
    
end
fprintf('\n')

end