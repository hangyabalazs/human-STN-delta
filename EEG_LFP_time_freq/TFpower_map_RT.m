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
    
    
    % Calculate and save correlation maps patient-by-patient
    TFmap_RTcorr_pats(sess2analyse,event,event,csd,bipol,baseline_win,chanmean)
    
    % Draw correlation maps patient-by-patient
    RTcorrmap_FIG_pats(sess2analyse,event,event,csd,bipol,baseline_win,plot_window,'bas',chanmean)
    
    % Draw correlation maps averaged across patients + perform statistics
    avgRTmap(sess2analyse,side,condi, event,event,freq_bands,plot_window,csd,'TrialNorm',baseline_win,chanmean)

    % Draw topoplots derived from correlation maps averaged across patients + perform statistics
    if strcmp(rectype, 'EEG') && strcmp(rectime, 'postop')
        topoplot_RTmap(sess2analyse,side,condi, event,event, fr_name,freq_bands,plot_window,...
            500,[-.3 .3],[-.1 .1],csd,bipol,'TrialNorm',false,true,baseline_win)
        
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
                
                
                corrmapfig(coefD,pvalD,normtype,times,timinx,f,patnm,condi,side,event,evty,act_chan,figdir);
            else
                fprintf('No channel data %s %s.\n',patnm, act_chan);
            end
            
            
        end
    elseif chanmean
        
        coefD = RTcorrMAP.([patnm '_' side]).chanmean.(event).(evty).(normtype).CoefMap;
        pvalD = RTcorrMAP.([patnm '_' side]).chanmean.(event).(evty).(normtype).PMap;
        
        
        corrmapfig(coefD,pvalD,normtype,times,timinx,f,patnm,condi,side,event,evty,'chanmean',figdir);
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

sess2analyse = sess2analyse(ismember({sess2analyse.side},side));  % safety check side and condition - plot only one side and one condition at once
sess2analyse = sess2analyse(ismember({sess2analyse.tag},condition));  % safety check side and condition - plot only one side and one condition at once


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
    side = sess2analyse(si).side;
    
    fprintf('%s, %s, %s...',patnm, condi, side)
    
    
    if ~chanmean
        chanlabs = fieldnames(RTcorrMAP.([patnm '_' side]));
        chnr = length(chanlabs)-1;
    else
        chanlabs = chanlabels;
        chnr = 1;
    end
    [coefdat, pdat] = deal(cell(1,chnr)); chansort = nan(1,chnr);
    
    for ic = 1:chnr
        act_chan = chanlabs{ic};
        coefdat{ic}= RTcorrMAP.([patnm '_' side]).(act_chan).(event).(evty).(normtype).CoefMap;
        pdat{ic}= RTcorrMAP.([patnm '_' side]).(act_chan).(event).(evty).(normtype).PMap;
        
        
    end
    
    coefdat2 = cat(3,coefdat{:});
    coefdat3 = coefdat2(f_ind,newtiminx,:);
    
    pdat2 = cat(3,pdat{:});
    maskdat = pdat2(f_ind,newtiminx,:)<0.05;
    
    
    
    allcoefs{si} = cat(3,coefdat3, nan(size(coefdat3,1), size(coefdat3,2), maxchnr - size(coefdat3,3)));
end

allcoef_cat = cat(4,allcoefs{:});

allcoef_stat = permute(allcoef_cat,[1 2 4 3]);
formula = 'mean(arg1,3);';



allcoef_stat_chanmean = nanmean(allcoef_stat,4);
[exactp_ersp,maskdat,~] = boostat_eeglab_J(allcoef_stat_chanmean,[],...
    alpha,iter,false,'fdr',formula,basinx);

fig = figure;
imagesc(newtimes,1:length(f_ind),mean(allcoef_stat_chanmean,3));
% crange = caxis;
crange = [-.1 .1]
hold on;
contour(newtimes,1:length(f_ind),maskdat,'Color','r')
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
