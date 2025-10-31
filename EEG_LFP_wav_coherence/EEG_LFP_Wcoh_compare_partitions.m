function EEG_LFP_Wcoh_compare_partitions(EventTypes,SubEventTypes,side,partition,patgr_nm,baseline_type)
% EEG_LFP_WCOH_COMPARE_PARTITIONS Compares partitioned trials of EEG-LFP coherence data.
% EEG_LFP_Wcoh_compare_partitions(...) gets and saves EEG-LFP wavelet coherence
% matrices for trials partitioned according to PARTITION, into SUBEVENTTYPES groups, aligned to
% EVENTTYPES, including SIDE sided tests, of patient group PATGR_NM. Data
% is normalized relative to baseline, defined by BASELINE_TYPE. 
% Permutation test with cluster-based correction is applied.
% Time-frequency difference maps are drawn and saved.
%
% See also: WCOH_AVG, STATCONDFIELDTRIP
%
% Johanna Petra Szabó, 10.2025
% Clinic for Neurosurgery and Neurointervention
% Semmelweis University, Budapest, Hungary
% szabo.johanna.petra@semmelweis.hu

global filesdir

% Parameters
plot_win = [-1 1];
baseline_win = [-1 -.5];
all_freqs = [1 80];
freqs = [1 4];

timewin_small = [-.5 0];
freq_small = [1.5 2.5];


if contains(patgr_nm, 'all')
    patients = 'allpatients';
else
%     patients = clinical_groups({patgr_nm},'intraop','stimoff','left');
    patients = clinical_groups({patgr_nm},'intraop','stimoff');
    patients = patients{1};
end
sess2analyse = getdata2analyse(filesdir, 'rectype','EEG',...
    'rectime','intraop','patients', patients, 'side',side, 'condition','stimoff');



% WCoherence patient-by-patient - average across event epochs (plot + save struct)
for ei = 1:length(EventTypes)

    event = EventTypes{ei}; fprintf('Event: %s,',event);

    for sei = 1:2
        subevent = SubEventTypes{ei,sei};
        downsamp = true;

        fprintf(' %s...\n',subevent);

%         wcoh_onebyone(sess2analyse,event,subevent,0,all_freqs,downsamp,NaN,plot_win,baseline_win,1)
        wcoh_onebyone(sess2analyse,event,subevent,0,all_freqs,downsamp,NaN,plot_win,[],1);

    end
end


% Compare partitions


for ei = 1:length(EventTypes)
    event = EventTypes{ei};
    subevs = SubEventTypes(ei,:);
    
    %% Compare wavelet coherence map (permutation test with cluster correction)
    compare_wcohmap_patientavg(event,subevs,plot_win,freqs,'all','all',0.05,'cluster',baseline_win,baseline_type,partition,patgr_nm,side)
    
    %% Compare selected time windows (paired T test)
    compare_partitions(event,subevs,freq_small,{[]},false,timewin_small,{},baseline_win,'all')
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
function compare_wcohmap_patientavg(event,subevs,plot_win,freqs,close2centr,...
    maxpower,alpha,mcorrect,baseline_win,baseline_type,partition,patgr_nm,side)

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
    'close2centr',close2centr,'maxpower',maxpower,'baseline_win',[],'plot_win',plot_win,'side',side);
    
    
basavg =  repmat(nanmean(basdat(:,basinx,:),2),[1 size(basdat,2) 1]);
bassd =  repmat(std(basdat(:,basinx,:),[],2,'omitnan'),[1 size(basdat,2) 1]);


statdat = {};
for ss = 1:2
    
    
    [sdat, ~,allpats{ss},ff] = wcoh_avg(event,'subevent',subevs{ss},...
        'region',{[]},'iscell',[],'freqs',freqs,'isfig',false,'downsamp',true,...
        'close2centr',close2centr,'maxpower',maxpower,'baseline_win',[],'plot_win',plot_win,'patients',patgr_nm);
    
    if strcmp(baseline_type,'indiv')
        basavg =  repmat(nanmean(sdat(:,basinx,:),2),[1 size(sdat,2) 1]);
        bassd =  repmat(std(sdat(:,basinx,:),[],2,'omitnan'),[1 size(sdat,2) 1]);
        
        sdat2 = (sdat-basavg)./ bassd;
        statdat{ss,1} = sdat2;
    else
        sdat_partss{ss} = sdat;
        patix = ismember(allpats_b,allpats{ss});
        sdat2 = (sdat-basavg(:,:,patix))./ bassd(:,:,patix);
        statdat{ss,1} = sdat2;
    end
    
%     
    
    
    sdat = [];
end

% if strcmp(baseline_type,'common')
%     sdat_all = mean( cat(4, sdat_partss{:}),4);
%     basavg =  repmat(nanmean(sdat_all(:,basinx,:,:),[2 4]),[1 size(sdat_all,2) 1]);
%     bassd =  repmat(std(sdat_all(:,basinx,:,:),[],[2 4],'omitnan'),[1 size(sdat_all,2) 1]);
%     
%     
%     statdat = {};
%     for ss = 1:2
%         
%         sdat2 = (sdat_partss{ss}-basavg)./ bassd;
%         statdat{ss,1} = sdat2;
%     end
% else
% end



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
    caxis([-1 1])
    title([subevs{sei}, ', n = ' num2str(length(allpats{sei}))])
end
subplot(3,1,3)
wcoh_fig(statscond{1},[],pltime,ff,[]); caxis([-2 2])
title([subevs{1} ' - ' subevs{2} ': T values'])

hold on; contour(pltime,log2(ff),pcond{1},'Color','red');

fnm = [partition(2:end) '_' side '_' event '_TIME' num2str(plot_win) '_FREQ' num2str(freqs) '_CLOSE2_' close2centr '_MAXPOW' maxpower '_ALPHA' num2str(alpha) '_' mcorrect '_BAS' num2str(baseline_win) '_' baseline_type '_' patgr_nm];
saveas(fig,fullfile(figdir_stat,[fnm '.jpg']))
saveas(fig,fullfile(figdir_stat,[fnm '.fig']))
% saveas(fig,fullfile(figdir_stat,[fnm '.pdf']))
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

