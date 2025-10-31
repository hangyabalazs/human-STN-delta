function TFpower_map_RT
% TFPOWER_MAP_RT  Correlation map between reaction time and wavelet power coefficients
%   TFPOWER_MAP_RT calculates trial-by-trial correlation (Pearson's correlation)
%   between reaction time and time-frequency wavelet power coefficients of each patient.
%   Significant correlation tested by TTEST on the level of individual patients
%   is marked by contouring on individual patient figures. Significant
%   change in correlation tested by permutation test wih fdr correction is
%   marked by contouring on patient-averaged figures.


% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

global filesdir

% Parameters
event = 'StimulusOn';
baseline_win = [-1 -.5];
plot_window = [-1 1];
bipol = false;
side = 'left';
condi = 'stimoff';

fr_name = 'high_delta';
freq_bands = [1 4];

for ri = 1:2
    switch ri;
        case 1; rectype = 'LFP'; rectime = 'intraop'; csd = false; chanmean = true;
        case 2; rectype = 'EEG'; rectime = 'postop'; csd = true; chanmean = false;
    end
    
    sess2analyse = getdata2analyse(filesdir, 'rectype',rectype,...
        'rectime',rectime,'patients', 'allpatients', 'side',side, 'condition',condi);
    
    %%
    % Calculate and save correlation maps patient-by-patient
    TFmap_RTcorr_pats(sess2analyse,event,event,csd,bipol,baseline_win,chanmean)
    
    %%
    % Draw correlation maps patient-by-patient (within-subject)
    RTcorrmap_FIG_pats(sess2analyse,event,event,csd,bipol,baseline_win,plot_window,'bas',chanmean)
    
    %%
    % Draw correlation maps averaged across patients + perform statistics
    avgRTmap(sess2analyse,side,condi, event,event,[1 80],plot_window,csd,'TrialNorm',baseline_win,chanmean)
%     avgRTmap(sess2analyse,side,condi, event,event,[1 80],plot_window,csd,'BASLims51_76',baseline_win,chanmean)
    
    %% Draw topoplots derived from correlation maps averaged across patients + perform statistics
    if strcmp(rectype, 'EEG') && strcmp(rectime, 'postop')
        topoplot_RTmap(sess2analyse,side,condi, event,event, fr_name,freq_bands,plot_window,...
            500,[-.3 .3],[-.1 .1],csd,bipol,'TrialNorm',false,true,baseline_win)
        
    end
    
    %% Calculate and draw between-subject correlation maps (median RT vs average TF)
    
    TFmap_RTcorr_betweenS(sess2analyse, plot_window, event,'rectype',rectype,...
        'chanmean',chanmean, 'csd',csd,'bipol',bipol,'baseline_win',baseline_win)
    
    TFmap_RTcorr_betweenS(sess2analyse, plot_window, event,'rectype',rectype,...
        'chanmean',chanmean, 'csd',csd,'bipol',bipol,'baseline_win',[])
    
    %%
    if strcmp(rectype, 'EEG') && strcmp(rectime, 'postop')
        topoplot_RTmap_betweenS(side,condi, event,event, fr_name,freq_bands,plot_window,...
            500,[-1 1],csd,bipol,'BASLims51_76')
        
        
%         topoplot_RTmap_betweenS(side,condi, event,event, fr_name,freq_bands,plot_window,...
%             500,[-.6 .1],csd,bipol,'TrialNorm')
%         
    end
end

end




%--------------------------------------------------------------------------
function TFmap_RTcorr_pats(sess2analyse,event,evty,csd,bipol,baseline_win,chanmean)

global rootdir figdir_pd

rectime = sess2analyse(1).rectime;
rectype = sess2analyse(1).rectype;

epoch_nm = [ event '_' evty '_EPOCHs.mat'];


figdir = fullfile(figdir_pd,[rectime '_' rectype],'behav_corr',...
    ['RT_corrmap_CSD' char(string(csd))]);

if bipol
    figdir = fullfile(figdir_pd,[rectime '_' rectype],'behav_corr',...
        ['RT_corrmap_CSD' char(string(csd)) '_bipol']);
end

if ~isfolder(figdir); mkdir(figdir); end;

if strcmp(rectype,'EEG')
    chanlocs = [];
    load(fullfile(rootdir,'postop_EEG_chanlocs.mat'));
    if bipol
        chanlabels = {'F4-F3'};
    else
        chanlabels = {chanlocs.labels};
    end
    if strcmp(rectime,'intraop')
        finx = strcmp('F4',chanlabels);
        chanlocs = chanlocs(finx);
        chanlabels = {'F4'};
    end
else
    chanlabels = arrayfun(@(x) ['Ch' num2str(x)] ,1:5,'UniformOutput',0);
end


%%

load(fullfile(rootdir,'freq_components.mat'));

epoch_win = [-2 2];
new_sr = 50;
times = epoch_win(1):1/new_sr:epoch_win(2);times = times(1:end-1);


% Baseline norm
baslims = dsearchn(times',baseline_win'); basinx = baslims(1):baslims(2);
basnm = ['BASLims' num2str(baslims(1)) '_' num2str(baslims(2))];

sessnr = length(sess2analyse);
for si = 1:sessnr
    
    if exist(fullfile(figdir,'RTcorrMAP.mat'))==2
        load(fullfile(figdir,'RTcorrMAP.mat'))
    else
        RTcorrMAP = struct;
    end
    
    curr_resdir = sess2analyse(si).folder;
    condi = sess2analyse(si).tag;
    patnm = sess2analyse(si).patient;
    side = sess2analyse(si).side;
    
    fprintf('%s, %s, %s...',patnm, condi, side)
    
    try
        load(fullfile(curr_resdir,'Evinxx.mat'));
    catch
        fprintf('No data %s %s %s\n',patnm, side,condi);
        continue;
    end
    
    evinx = Evinxx.(event).(evty).TE_index;
    
    % Reaction Time
    
    [RT,Go_RT,FAlarm_RT,Hit,NoStopTrials StopTrials,perf,stopperf] = ...
        calc_RT(sess2analyse(si),['TrialEvents_' condi '.mat'],false);
    RT2 = RT{1}(evinx);
    
    
    %% Delta map
    
    % Load norm. maps
    %     if bipol
    %         resdir = fullfile(curr_resdir, 'EventAVGs_bipol');
    %     else
    %         resdir = fullfile(curr_resdir, ['EventAVGs_CSD' char(string(csd))]);
    %     end
    %
    %     load(fullfile(resdir,epoch_nm))
    %
    %     chanlabs = fieldnames(EventEPOCHs);
    
    
    if csd
        tfdirnm = 'TFpows_blocks_CSD';
    elseif bipol
        tfdirnm = 'TFpows_blocks_bipol';
    else
        tfdirnm = 'TFpows_blocks';
    end
    
    
    
    
    
    
    % Loopp over channels
    chnr = length(chanlabels);
    evpow_allc = cell(chnr,1);
    
    for ic = 1:chnr
        
        
        act_chan = chanlabels{ic};
        fprintf('%s...',act_chan)
        %                 tf_1chan = EventEPOCHs.(act_chan).TrialNorm;
        
        
        % Load TF blocks
        [epoch_pow, f] = load_TFblocks(curr_resdir,act_chan,'Pow',fullfile(sess2analyse(si).folder,tfdirnm),4);
        if isempty(epoch_pow)
            fprintf('No channel data %s %s.\n',patnm, act_chan);
            continue
        end
        
        if ~strcmp(event,'StopSignal') && ismember(evty,{'FailedStopTrial','SuccesfulStopTrial'});
            
            [~,evinx] = StimOn_stoppart_evinx(Evinxx,event,evty);
        else
            evinx = Evinxx.(event).(evty).epoch_index;
        end
        
        if isempty(evinx)
            fprintf('No %s trial.\n',evty);
            continue
        end
        
        try
            evpow = epoch_pow(:,:,evinx);
        catch
            fprintf('evinx mismatch %s \n', curr_resdir)
        end
        
        if bipol
            act_chan(strfind(act_chan, '-')) = '_';
        end
        
        if ~chanmean
            % Full-epoch length single-trial norm. (TrialNorm)
            M = mean(evpow,2);
            SD = std(evpow,[],2);
            Mrep = repmat(M,[1 size(evpow,2)]);
            SDrep = repmat(SD,[1 size(evpow,2)]);
            tf_1chan = (evpow - Mrep)./SDrep;
            
            if size( tf_1chan ,3)~= length(evinx)
                fprintf('Trnr problem\n');
                keyboard;
            end
            
            RTcorrMAP = RTmap_save(RT2,tf_1chan,'TrialNorm',patnm,side,event,evty,act_chan,figdir,RTcorrMAP);
            
            
            
            B2 = repmat(nanmean(evpow(:,basinx,:),2),[1 200]);
            Bsd2 = repmat(std(evpow(:,basinx,:),[],2,'omitnan'),[1 200]);
            
            tf_1chan = (evpow- B2)./ Bsd2;
            
            RTcorrMAP = RTmap_save(RT2,tf_1chan,basnm,patnm,side,event,evty,act_chan,figdir,RTcorrMAP);
        end
        evpow_allc{ic} = evpow;
        
    end
    
    if chanmean
        
        evpow = mean( cat(4,evpow_allc{:}) ,4);
        
        
        M = mean(evpow,2);
        SD = std(evpow,[],2);
        Mrep = repmat(M,[1 size(evpow,2)]);
        SDrep = repmat(SD,[1 size(evpow,2)]);
        tf_1chan = (evpow - Mrep)./SDrep;
        
        RTcorrMAP = RTmap_save(RT2,tf_1chan,'TrialNorm',patnm,side,event,evty,'chanmean',figdir,RTcorrMAP);
        
        
        B2 = repmat(nanmean(evpow(:,basinx,:),2),[1 200]);
        Bsd2 = repmat(std(evpow(:,basinx,:),[],2,'omitnan'),[1 200]);
        
        tf_1chan = (evpow- B2)./ Bsd2;
        
        RTcorrMAP = RTmap_save(RT2,tf_1chan,basnm,patnm,side,event,evty,'chanmean',figdir,RTcorrMAP);
        
    end
    fprintf('\n')
    
end
end


%--------------------------------------------------------------------------
function RTcorrMAP = RTmap_save(RT2,tf_1chan,normtype,patnm,side,event,evty,act_chan,figdir,RTcorrMAP)
rtnn = ~isnan(RT2);

tf = permute(tf_1chan(:,:,rtnn),[3 1 2]);

if size(RT2,1)==1; RT2 = RT2'; end;
rtmat = repmat(RT2(rtnn), [1 size(tf_1chan,[1 2])] );

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
function RTcorrmap_FIG_pats(sess2analyse,event,evty,csd,bipol,baseline_win,plot_window,norm,chanmean)

global rootdir figdir_pd

rectime = sess2analyse(1).rectime;
rectype = sess2analyse(1).rectype;



figdir = fullfile(figdir_pd,[rectime '_' rectype],'behav_corr',...
    ['RT_corrmap_CSD' char(string(csd))]);

if bipol
    figdir = fullfile(figdir_pd,[rectime '_' rectype],'behav_corr',...
        ['RT_corrmap_CSD' char(string(csd)) '_bipol']);
end


if strcmp(rectype,'EEG')
    chanlocs = [];
    load(fullfile(rootdir,'postop_EEG_chanlocs.mat'));
    if bipol
        chanlabels = {'F4-F3'};
    else
        chanlabels = {chanlocs.labels};
    end
    if strcmp(rectime,'intraop')
        finx = strcmp('F4',chanlabels);
        chanlocs = chanlocs(finx);
        chanlabels = {'F4'};
    end
else
    chanlabels = arrayfun(@(x) ['Ch' num2str(x)] ,1:5,'UniformOutput',0);
end


%%

load(fullfile(rootdir,'freq_components.mat'));

epoch_win = [-2 2];
new_sr = 50;
times = epoch_win(1):1/new_sr:epoch_win(2);times = times(1:end-1);
newtimes = plot_window(1):1/new_sr:plot_window(2);newtimes = newtimes(1:end-1);
timinx = dsearchn(times',newtimes');

if contains(norm,'bas');
    baslims = dsearchn(times',baseline_win');
    normtype = ['BASLims' num2str(baslims(1)) '_' num2str(baslims(2))];
else
    normtype = 'TrialNorm';
end

load(fullfile(figdir,'RTcorrMAP.mat'))

sessnr = length(sess2analyse);
for si = 1:sessnr
    
    
    condi = sess2analyse(si).tag;
    patnm = sess2analyse(si).patient;
    side = sess2analyse(si).side;
    
    fprintf('%s, %s, %s...',patnm, condi, side)
    
    
    if ~chanmean
        % Loopp over channels
        chnr = length(chanlabels);
        
        for ic = 1:chnr
            
            act_chan = chanlabels{ic};
            fprintf('%s...',act_chan)
            
            if isfield(RTcorrMAP.([patnm '_' side]),act_chan)
                coefD = RTcorrMAP.([patnm '_' side]).(act_chan).(event).(evty).(normtype).CoefMap;
                pvalD = RTcorrMAP.([patnm '_' side]).(act_chan).(event).(evty).(normtype).PMap;
                
                
                corrmapfig(coefD,pvalD,normtype,times,timinx,f,patnm,condi,side,event,evty,act_chan,figdir,'none');
            else
                fprintf('No channel data %s %s.\n',patnm, act_chan);
            end
            
            
        end
    elseif chanmean
        
        coefD = RTcorrMAP.([patnm '_' side]).chanmean.(event).(evty).(normtype).CoefMap;
        pvalD = RTcorrMAP.([patnm '_' side]).chanmean.(event).(evty).(normtype).PMap;
        
        
        corrmapfig(coefD,pvalD,normtype,times,timinx,f,patnm,condi,side,event,evty,'chanmean',figdir,'none',rectime,rectype);
    end
    fprintf('\n')
    
end
end


%--------------------------------------------------------------------------
function topoplot_RTmap(sess2analyse,side,condition, event,evty, fr_name,freq_bands,plot_win,...
    topobin,pat_cLim,avg_cLim,csd,bipol,normtype,patfig,avgfig,baseline_win)

global figdir_pd rootdir

rectime = sess2analyse(1).rectime;
rectype = sess2analyse(1).rectype;

sess2analyse = sess2analyse(ismember({sess2analyse.side},side));  % safety check side and condition - plot only one side and one condition at once
sess2analyse = sess2analyse(ismember({sess2analyse.tag},condition));  % safety check side and condition - plot only one side and one condition at once


load(fullfile(rootdir,'freq_components.mat'));
f_ind = intersect(find(f>=freq_bands(1)),find(f<=freq_bands(2)));

alpha = 0.05;
iter = 1000;


new_sr = 50;
newtimes = plot_win(1):1/new_sr:plot_win(2);newtimes = newtimes(1:end-1);
baslims = dsearchn(newtimes',baseline_win'); basinx = baslims(1):baslims(2);


if strcmp(rectype,'EEG')
    chanlocs = [];
    load(fullfile(rootdir,'postop_EEG_chanlocs.mat'));
    if bipol
        chanlabels = {'F4-F3'};
    else
        chanlabels = {chanlocs.labels};
    end
    if strcmp(rectime,'intraop')
        finx = strcmp('F4',chanlabels);
        chanlocs = chanlocs(finx);
        chanlabels = {'F4'};
    end
else
    chanlabels = arrayfun(@(x) ['Ch' num2str(x)] ,1:5,'UniformOutput',0);
end

maxchnr = length(chanlabels);

figdir = fullfile(figdir_pd,[rectime '_' rectype],'behav_corr',...
    ['RT_corrmap_CSD' char(string(csd))]);

if bipol
    figdir = fullfile(figdir_pd,[rectime '_' rectype],'behav_corr',...
        ['RT_corrmap_CSD' char(string(csd)) '_bipol']);
end
if contains(maptype,'between')
    figdir = [figdir '_betweenSubj'];
end



epoch_win = [-2 2];
new_sr = 50;
times = epoch_win(1):1/new_sr:epoch_win(2);times = times(1:end-1);
newtimes = plot_win(1):1/new_sr:plot_win(2);newtimes = newtimes(1:end-1);
newtiminx = dsearchn(times',newtimes');


load(fullfile(figdir,'RTcorrMAP.mat'));


allcoefs = cell(1,length(sess2analyse));
for si = 1:length(sess2analyse)
    condi = sess2analyse(si).tag;
    patnm = sess2analyse(si).patient;
    side = sess2analyse(si).side;
    
    fprintf('%s, %s, %s...',patnm, condi, side)
    
    
    
    chanlabs = fieldnames(RTcorrMAP.([patnm '_' side]));
    chnr = length(chanlabs);
    [coefdat, pdat] = deal(cell(1,chnr)); chansort = nan(1,chnr);
    
    for ic = 1:chnr
        act_chan = chanlabs{ic};
        coefdat{ic}= RTcorrMAP.([patnm '_' side]).(act_chan).(event).(evty).(normtype).CoefMap;
        pdat{ic}= RTcorrMAP.([patnm '_' side]).(act_chan).(event).(evty).(normtype).PMap;
        
        
        
        if strcmp(rectype,'EEG') && strcmp(rectime,'postop')
            chansort(ic) = find( ismember({chanlocs.labels}, act_chan) );
        end
    end
    
    coefdat2 = cat(3,coefdat{:});
    coefdat3 = coefdat2(f_ind,newtiminx,:);
    
    pdat2 = cat(3,pdat{:});
    maskdat = pdat2(f_ind,newtiminx,:)<0.05;
    
    
    if patfig
        figdir2 = fullfile(figdir,[patnm '_' side],[event '_' evty]); if ~isfolder(figdir2);mkdir(figdir2); end;
        
        if strcmp(rectype,'EEG') && strcmp(rectime,'postop')
            topoplot_fig(coefdat3,chanlocs(chansort),new_sr,plot_win,topobin,pat_cLim,figdir2,...
                maskdat,fr_name,[patnm '_' side '_' condi])
            
        end
        
        
    end
    
    allcoefs{si} = cat(3,coefdat3, nan(size(coefdat3,1), size(coefdat3,2), maxchnr - size(coefdat3,3)));
end

allcoef_cat = cat(4,allcoefs{:});
allcoef_avg = nanmean( allcoef_cat, 4 );

allcoef_stat = permute(allcoef_cat,[1 2 4 3]);
formula = 'mean(arg1,3);';



if avgfig
    if strcmp(rectype,'EEG') && strcmp(rectime,'postop')
        
        for ci = 1:chnr
            [exactp_ersp,maskersp{ci},~] = boostat_eeglab_J(allcoef_stat(:,:,:,ci),...
                [],0.01,iter,false,'fdr',formula,basinx);
        end
        
        maskdat = cat(3,maskersp{:});
        topoplot_fig(allcoef_avg,chanlocs(chansort),new_sr,plot_win,topobin,avg_cLim,figdir,...
            maskdat,fr_name,['AVG_' side '_' condi '_' event '_' evty])
    else
        
        allcoef_stat_chanmean = nanmean(allcoef_stat,4);
        [exactp_ersp,maskdat,~] = boostat_eeglab_J(allcoef_stat_chanmean,[],...
            alpha,iter,false,'fdr',formula,basinx);
        
        fig = figure;
        imagesc(newtimes,1:length(f_ind),mean(allcoef_stat_chanmean,3));
        crange = caxis;
        hold on;
        contour(newtimes,1:length(f_ind),maskdat,'Color',[.5 .5 .5])
        caxis(crange)
        yti = round(linspace(1,length(f_ind), 4)) ;
        yticks(yti); yticklabels( arrayfun(@num2str,f(f_ind(yti)),'UniformOutput',0) );
        
        colormap(bone);
        title({['Patietn AVG, channelmean'],[side ', ' condition]});
        xlabel('Time (s)'); ylabel('Frequency (Hz)');
        c=colorbar; ylabel(c,'R','Rotation',270);
        fnm = fullfile(figdir,['AVG_' side '_' condi '_' event '_' evty 'FR' num2str(freq_bands)]);
        saveas(fig,[fnm '.jpg'])
        saveas(fig,[fnm '.fig'])
        saveas(fig,[fnm '.pdf'])
        close(fig);
    end
end

end




%--------------------------------------------------------------------------
function avgRTmap(sess2analyse,side,condition, event,evty,freq_bands,plot_win,csd,normtype,baseline_win,chanmean)

global figdir_pd rootdir

rectime = sess2analyse(1).rectime;
rectype = sess2analyse(1).rectype;

if ~contains(side, 'both')
    sess2analyse = sess2analyse(ismember({sess2analyse.side},side));  % safety check side and condition - plot only one side and one condition at once
end
if ~contains(condition, 'both')
    sess2analyse = sess2analyse(ismember({sess2analyse.tag},condition));  % safety check side and condition - plot only one side and one condition at once
end

load(fullfile(rootdir,'freq_components.mat'));
f_ind = intersect(find(f>=freq_bands(1)),find(f<=freq_bands(2)));

alpha = 0.05;
iter = 1000;


new_sr = 50;
newtimes = plot_win(1):1/new_sr:plot_win(2);newtimes = newtimes(1:end-1);
baslims = dsearchn(newtimes',baseline_win'); basinx = baslims(1):baslims(2);


if strcmp(rectype,'EEG') && strcmp(rectime,'intraop')
    chanlabels = {'F4'};
else
    if ~chanmean
        chanlabels = arrayfun(@(x) ['Ch' num2str(x)] ,1:5,'UniformOutput',0);
    else
        chanlabels = {'chanmean'};
    end
end

maxchnr = length(chanlabels);

figdir = fullfile(figdir_pd,[rectime '_' rectype],'behav_corr',...
    ['RT_corrmap_CSD' char(string(csd))]);


epoch_win = [-2 2];
new_sr = 50;
times = epoch_win(1):1/new_sr:epoch_win(2);times = times(1:end-1);
newtimes = plot_win(1):1/new_sr:plot_win(2);newtimes = newtimes(1:end-1);
newtiminx = dsearchn(times',newtimes');


load(fullfile(figdir,'RTcorrMAP.mat'));


allcoefs = cell(1,length(sess2analyse));
for si = 1:length(sess2analyse)
    condi = sess2analyse(si).tag;
    patnm = sess2analyse(si).patient;
    pside = sess2analyse(si).side;
    
    fprintf('%s, %s, %s...',patnm, condi, pside)
    
    
    if ~chanmean
        chanlabs = fieldnames(RTcorrMAP.([patnm '_' pside]));
        chnr = length(chanlabs)-1;
    else
        chanlabs = chanlabels;
        chnr = 1;
    end
    [coefdat, pdat] = deal(cell(1,chnr)); chansort = nan(1,chnr);
    
    for ic = 1:chnr
        act_chan = chanlabs{ic};
        coefdat{ic}= RTcorrMAP.([patnm '_' pside]).(act_chan).(event).(evty).(normtype).CoefMap;
        pdat{ic}= RTcorrMAP.([patnm '_' pside]).(act_chan).(event).(evty).(normtype).PMap;
        
        
    end
    
    coefdat2 = cat(3,coefdat{:});
    coefdat3 = coefdat2(f_ind,newtiminx,:);
    
    pdat2 = cat(3,pdat{:});
    maskdat = pdat2(f_ind,newtiminx,:)<0.05;
    
    
    
    allcoefs{si} = cat(3,coefdat3, nan(size(coefdat3,1), size(coefdat3,2), maxchnr - size(coefdat3,3)));
end

allcoef_cat = cat(4,allcoefs{:});

allcoef2 = permute(allcoef_cat,[1 2 4 3]);
% formula = 'mean(arg1,3);';



allcoef_stat = nanmean(allcoef2,4); % channel mean


if contains(lower(normtype), 'trial')
    allcoef_statN = allcoef_stat;
    allcoef_statMean = mean(allcoef_stat,3);
else
    %     bas_avgx = repmat(mean(allcoef_stat(:,basinx,:),2), [1 size(allcoef_stat,2) 1]);
    %     bas_stdx = repmat(std(allcoef_stat(:,basinx,:),[], 2,'omitnan'), [1 size(allcoef_stat,2 ) 1]);
    %     allcoef_statN = (allcoef_stat - bas_avgx)./bas_stdx;
    bas_avg = repmat(mean(allcoef_stat(:,basinx,:),[2 3]), [1 size(allcoef_stat,2) 1]);
    bas_std = repmat(std(allcoef_stat(:,basinx,:),[], [2 3],'omitnan'), [1 size(allcoef_stat,2) 1]);
    allcoef_statMean = (mean(allcoef_stat,3) - bas_avg)./bas_std;
end
formula = 'mean(arg1,3);';
[exactp_ersp,maskdat,~] = boostat_eeglab_J(allcoef_stat,[],alpha,iter,false,'fdr',formula,basinx);
  
fig = figure;
imagesc(newtimes,1:length(f_ind),allcoef_statMean);

% crange = caxis;
crange = [-.1 .1];
hold on;
zmap = norminv(exactp_ersp);
bootstatFDR_clustercorr(zmap,maskdat,f(f_ind),rectime,rectype,act_chan,newtimes,1:length(f_ind))

% contour(newtimes,1:length(f_ind),maskdat,'Color','r')
caxis(crange)
yti = round(linspace(1,length(f_ind), 4)) ;
yticks(yti); yticklabels( arrayfun(@num2str,f(f_ind(yti)),'UniformOutput',0) );

colormap(bone);
title({['Patietn AVG, channelmean'],[side ', ' condition]});
xlabel('Time (s)'); ylabel('Frequency (Hz)');
c=colorbar; ylabel(c,'R','Rotation',270);
fnm = fullfile(figdir,['AVG_' side '_' condi '_' event '_' evty 'FR' num2str(freq_bands) '_' normtype]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);


end



%---------------------------------------------------------------------------
function TFmap_RTcorr_betweenS(sess2analyse, plot_win,event,varargin)


prs = inputParser;
addRequired(prs,'sess2analyse',@isstruct);
addRequired(prs,'plot_win',@isvector);
addRequired(prs,'event',@ischar);
addParameter(prs,'rectype','',@ischar);
addParameter(prs,'epoch_win',[-2 2],@isvector);
addParameter(prs,'new_sr',50,@isnumeric);
addParameter(prs,'chanmean',1,@(x) isnumeric(x)|islogical(x));
addParameter(prs,'csd',false,@islogical);
addParameter(prs,'bipol',false,@islogical);
addParameter(prs,'baseline_win',[],@(x) isvector(x)|isnumeric(x));
% addParameter(prs,'pow_normtype','TRN_AVG',@ischar);

parse(prs,sess2analyse,plot_win,event,varargin{:});
p = prs.Results;

global rootdir figdir_pd


% Recordint time
rectime = sess2analyse(1).rectime;
if ~strcmp(rectime,'postop')
    p.csd = false; p.bipol = false;
end

condi = sess2analyse(1).tag;

% Result directory
figdir = fullfile(figdir_pd,[rectime '_' p.rectype],'behav_corr',...
    ['RT_corrmap_CSD' char(string(p.csd))]);
if p.bipol;  figdir = [figdir '_bipol']; end

figdir = [figdir '_betweenSubj'];
if ~isdir(figdir); mkdir(figdir); end;


% Load freq components
load(fullfile(rootdir,'freq_components.mat'))

% Time

if ~isequal(p.plot_win,p.epoch_win)
    times = p.epoch_win(1):1/p.new_sr:p.epoch_win(2);
    new_times = p.plot_win(1):1/p.new_sr:p.plot_win(2);
    newtinx = dsearchn(times',new_times'); newtinx = newtinx(1:end-1);
else
    newtinx = 1:200;
end

if ~isempty(p.baseline_win)
    
    baslims = dsearchn(times',p.baseline_win');
%     basinx = baslims(1):baslims(2);
    basnm = ['BASLims' num2str(baslims(1)) '_' num2str(baslims(2))]; normtype = basnm;
    
    if newtinx(1)>=baslims(1)
        newtinx = newtinx-baslims(1)+1; % bas-corrected data is saved from start of baseline until end of original data epoch
    else
        fprintf('Start of time period before baseline start\n');
        return;
    end
    timerang = [max( [ p.baseline_win(1) p.plot_win(1) ]) p.plot_win(2)];
else
    basnm = ['TRN_AVG'];
    normtype = 'TrialNorm';
    timerang = p.plot_win;
end

avgTF = {};
sesnr = length(sess2analyse);
for snr = 1:sesnr
    curr_resdir = sess2analyse(snr).folder;
    side = sess2analyse(snr).side;
    patnm = sess2analyse(snr).patient;
    
    
    fprintf('%s %s %s...\n', patnm,side);
    
    % Load event averages for patients
    if p.bipol
        resdir = fullfile(curr_resdir, 'EventAVGs_bipol');
    else
        resdir = fullfile(curr_resdir, ['EventAVGs_CSD' char(string(p.csd))]);
    end
    
    
    
    
    % Channels
    if strcmp(p.rectype,'LFP') && p.chanmean==1
        channels = {'chanmean'};
    elseif strcmp(p.rectype,'EEG') && strcmp(rectime,'intraop')
        channels = {'F4'}; p.chanmean = 0;
    elseif strcmp(p.rectype,'EEG') && strcmp(rectime,'postop')
        if ~p.bipol
            load(fullfile(rootdir,'postop_EEG_chanlocs.mat'))
            channels = {chanlocs.labels};
        else
            channels = {'F4-F3'};
        end
    end
    
    if isempty(avgTF)
        avgTF = cell(length(sess2analyse),length(channels));
    end
    
    if strcmp(p.rectype,'LFP') && p.chanmean==1
        
        evavg_nm = [ event '_' event '_chan_AVGs.mat'];
        evSTAT_nm = [ event '_' event '_chan_STATs.mat'];
    else
        evavg_nm = [ event '_' event '_AVGs.mat'];
        evSTAT_nm = [ event '_' event '_STATs.mat'];
    end
    
    
    try
        load(fullfile(resdir,evavg_nm))
        load(fullfile(resdir,evSTAT_nm))
    catch
        fprintf('No files %s %s %s %s %s \n', patnm, side, event)
        continue;
    end
    
    
    
    if strcmp(p.rectype,'LFP') && p.chanmean~=1
        channels = fieldnames(EventAVGs);
    end
    
    % Loop over channels
    for ch = 1:length(channels)
        act_chan = channels{ch};
        if p.bipol
            act_chan(strfind(act_chan, '-')) = '_';
        end
        
        
        
        
        tf = EventAVGs.(act_chan).(basnm);
        
        avgTF{snr,ch} = tf;
        
    end
    
    
end


% Median RT
rt_figdir = fullfile(figdir_pd, 'Behav','RTmedian');
load(fullfile(rt_figdir,['All_RTs_all.mat']));

medRT = nan(length(sess2analyse),1 );
for snr = 1:sesnr
    
    patnm = sess2analyse(snr).patient;
    tag = sess2analyse(snr).tag;
    side = sess2analyse(snr).side;
    
    rt_patinx =  find(strcmp(patnm,{All_RTs.patient}));
    rt_inx = find(strcmp(side,{All_RTs(rt_patinx).side}));
    
    medRT(snr) = nanmedian(rmoutliers(All_RTs(rt_patinx(rt_inx)).([rectime '_' tag])));
    
end


if length( unique({sess2analyse.side}) )==2; sidnm = 'bothside'; else; sidnm = side; end;

% Correlate

if exist(fullfile(figdir,'RTcorrMAP.mat'))==2
    load(fullfile(figdir,'RTcorrMAP.mat'))
else
    RTcorrMAP = struct;
end


chnr = size(avgTF,2);
patnm = 'Between_patients';
for ch = 1:chnr
    
    act_chan = channels{ch};
    
    
    medRT2 = medRT; medRT2(cellfun(@isempty,avgTF(:,ch))) = [];
    tf_1chan = cat(3,avgTF{:,ch});
    RTcorrMAP = RTmap_save(medRT2, tf_1chan,normtype,patnm,sidnm,event,event,act_chan,figdir,RTcorrMAP);
    
    coefD = RTcorrMAP.([patnm '_' sidnm]).(act_chan).(event).(event).(normtype).CoefMap ;
    pvalD =  RTcorrMAP.([patnm '_' sidnm]).(act_chan).(event).(event).(normtype).PMap;
    
%     if ~isempty(basline_win)
%     coefBa = mean (mean( coefD(:,basinx,:),2 ) ,3);
%     coefBs = std (mean( coefD(:,basinx,:),2 ) ,[],3);
%     end
    
    corrmapfig(coefD,pvalD,normtype,times,newtinx,f,patnm,condi,sidnm,event,event,...
        act_chan,figdir,'none');
    
    
end


end




%--------------------------------------------------------------------------
function topoplot_RTmap_betweenS(side,condition, event,evty, fr_name,freq_bands,plot_win,...
    topobin,pat_cLim,csd,bipol,normtype)

global figdir_pd rootdir

load(fullfile(rootdir,'freq_components.mat'));
f_ind = intersect(find(f>=freq_bands(1)),find(f<=freq_bands(2)));

alpha = 0.05;
iter = 1000;


new_sr = 50;
newtimes = plot_win(1):1/new_sr:plot_win(2);newtimes = newtimes(1:end-1);


rectype = 'EEG'; rectime = 'postop';
chanlocs = [];
load(fullfile(rootdir,'postop_EEG_chanlocs.mat'));

chanlabels = {chanlocs.labels};



figdir = fullfile(figdir_pd,[rectime '_' rectype],'behav_corr',...
    ['RT_corrmap_CSD' char(string(csd))]);


if bipol
    figdir = fullfile(figdir_pd,[rectime '_' rectype],'behav_corr',...
        ['RT_corrmap_CSD' char(string(csd)) '_bipol']);
end

figdir = [figdir '_betweenSubj'];




epoch_win = [-2 2];
new_sr = 50;
times = epoch_win(1):1/new_sr:epoch_win(2);times = times(1:end-1);
newtimes = plot_win(1):1/new_sr:plot_win(2);newtimes = newtimes(1:end-1);
newtiminx = dsearchn(times',newtimes');


load(fullfile(figdir,'RTcorrMAP.mat'));


chanlabs = fieldnames(RTcorrMAP.(['Between_patients_' side]));
chnr = length(chanlabs);
[coefdat, pdat] = deal(cell(1,chnr)); chansort = nan(1,chnr);

for ic = 1:chnr
    act_chan = chanlabs{ic};
    coefdat{ic}= RTcorrMAP.(['Between_patients_' side]).(act_chan).(event).(evty).(normtype).CoefMap;
    pdat{ic}= RTcorrMAP.(['Between_patients_' side]).(act_chan).(event).(evty).(normtype).PMap;
    
    
    
  
        chansort(ic) = find( ismember({chanlocs.labels}, act_chan) );
    
end

coefdat2 = cat(3,coefdat{:});
coefdat3 = coefdat2(f_ind,newtiminx,:);

pdat2 = cat(3,pdat{:});
newalpha = fdr(pdat2,0.05);
maskdat = pdat2(f_ind,newtiminx,:)<=newalpha;


figdir2 = fullfile(figdir,[event '_' evty]); if ~isfolder(figdir2);mkdir(figdir2); end;

topoplot_fig(coefdat3,chanlocs(chansort),new_sr,plot_win,topobin,pat_cLim,figdir2,...
    maskdat,fr_name,['Topoplot_' side '_' condition '_' normtype])





end
