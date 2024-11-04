function updrs_behav_corr(conditions)
%UPDRS_BEHAV_CORR calculates correlation between UPDRS scores and behav. parameters
%   UPDRS_BEHAV_CORR calculates Pearson's correalation between preop and postop (DBS-on)
%   UPDRS scores (tested both with and without medication) and reaction time (RT), stop signal delay at 50%
%   probability of stopping (SSDp0.5), performance and stop performance.
%   Creates scatter plots and overlays regression line (robust regression) with confidence interval.
%
% Input parameters:
%       CONDITIONS      Nx2 cell array of task conditions to compare
%               -first column corresponds to recording time, second column to
%               DBS stimulation condition
%               - UPDRS scores are only calculated for preop and postop
%               (DBS-on) conditions
%               -ex.:{'preop','stimoff';'postop','stimon'}; 

% See also POLYPREDCICALL_MOD, CALC_RT, GETDATA2ANALYSE

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

global rootdir figdir_pd

if exist(fullfile(rootdir,'UPDRS_scores.mat'))~=2
    % Create UPRS tables (manually from excel table)
    fprintf('Create a cell array manually with patient codes with variable name: patients !')
    keyboard;
    fprintf('Create following variables (column vectors) manually with clinical parameters: \n \n preop_medOFF, preop_medON, postop_medOFF, postop_medON,year_of_onset,disease_duration,postop_motor_impr,LEDD_reduction !')
    keyboard;
    updrs_tab = table(preop_medOFF,preop_medON,postop_medOFF,postop_medON,...
        year_of_onset,disease_duration,postop_motor_impr,LEDD_reduction,'RowNames',patients);
    save(fullfile(rootdir,'UPDRS_scores.mat'),'updrs_tab')
end


figdir = fullfile(figdir_pd,'Behav','UPDRS_corrs'); if ~isfolder(figdir); mkdir(figdir); end;
condnr = size(conditions,1);
for coci = 1:condnr; % loop over conditions
    rectime = conditions{coci,1};
    condi = conditions{coci,2};
    
    
    % RT
    
    for si = 1:2
        switch si; case 1; side = 'left'; case 2; side = 'right'; end;
        rt_updrs_corr(rectime, side,condi,figdir)
    end
    
    
    
    % SSD0.5
    
    for si = 1:2
        switch si; case 1; side = 'left'; case 2; side = 'right'; end;
        ssrt_updrs_corr(rectime, side,condi,figdir,'ssd05')
        
    end
    
    
    % Performance
    
    for si = 1:2
        switch si; case 1; side = 'left'; case 2; side = 'right'; end;
        perf_updrs_corr(rectime, side,condi,figdir)
    end
    
end



end





%--------------------------------------------------------------------------
function rt_updrs_corr(rectime, side,condi,figdir)

global rootdir

load(fullfile(rootdir,'UPDRS_scores.mat'))

SessList = getdata2analyse(rootdir, 'rectype','BEHAV','rectime',rectime, ...
    'patients','allpatients', 'side',side, 'condition',condi);


[RT,Go_RT,FAlarm_RT,Hit,NoStopTrials StopTrials,perf,stopperf] = calc_RT(SessList,['TrialEvents_nosync_' condi '.mat'],false);
RT2 = cellfun(@(x) nanmedian(rmoutliers(x)),RT,'UniformOutput',false);

pats = {SessList.patient};
tabpats = updrs_tab.Properties.RowNames;
[~,ai,bi] = intersect(pats,tabpats);


fig = figure;

%  MED OFF
subplot(121)
upd = updrs_tab{bi,[rectime '_medOFF']};
rts = cell2mat(RT2(ai))';
% [p, R] = polypredcicall_mod(upd,rts,0.95,'robust',0.1);
[p, R] = polypredcicall_mod(upd,rts,0.95,'',0.1)
hold on; scatter(upd, rts,[],'k','filled');
xlabel([rectime ' UPDRS med OFF']); ylabel('RT (s)');
title([rectime ' med OFF vs median RT'])


%  MED ON
subplot(122)
upd = updrs_tab{bi,[rectime '_medON']};
rts = cell2mat(RT2(ai))';
% [p, R] = polypredcicall_mod(upd,rts,0.95,'robust',0.1);
[p, R] = polypredcicall_mod(upd,rts,0.95,'',0.1)
hold on; scatter(upd, rts,[],'k','filled');
xlabel([rectime ' UPDRS med ON']); ylabel('RT (s)');
title([rectime ' med ON vs median RT'])
suptitle([rectime ' ' side ' ' condi])

fnm = fullfile(figdir,[upper(rectime) '_medRT_' side '_' condi]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);

end


%--------------------------------------------------------------------------
function ssrt_updrs_corr(rectime, side,condi,figdir,param)

global rootdir figdir_pd

load(fullfile(rootdir,'UPDRS_scores.mat'))


ssrt_dir = fullfile(figdir_pd,'Behav','SSRT');
load(fullfile(ssrt_dir,[rectime '_' condi],'SSRT_results.mat'));

pats = fieldnames(SSRT_results);
tabpats = updrs_tab.Properties.RowNames;
[~,ai,bi] = intersect(pats,tabpats);

ss = nan(1,length(pats));
for k = 1:length(pats)
    
    sess = fieldnames(SSRT_results.(pats{k}));
    sinx = contains(sess,side(1));
    if any(sinx)
        ss(k) = SSRT_results.(pats{k}).(sess{sinx}).(param);
    end
end
ss = ss/1000;

% Remove wrong SSDp0.5 values
ss(ss<=0) = NaN;

fig = figure;

%  MED OFF
subplot(121)
upd = updrs_tab{bi,[rectime '_medOFF']};
ss2 = ss(ai);
nn = ~isnan(ss2);
ss2 = ss2(nn);
upd = upd(nn);
% [p, R] = polypredcicall_mod(upd,ss2,0.95,'robust',0.1);
[p, R] = polypredcicall_mod(upd,ss2,0.95,'',0.1)
hold on; scatter(upd, ss2,[],'k','filled');
xlabel([rectime ' UPDRS med OFF']); ylabel([param ' (s)']);
title([rectime ' med OFF vs ' param])


%  MED ON
subplot(122)
upd = updrs_tab{bi,[rectime '_medON']};
ss2 = ss(ai);
nn = ~isnan(ss2);
ss2 = ss2(nn);
upd = upd(nn);
% [p, R] = polypredcicall_mod(upd,ss2,0.95,'robust',0.1);
[p, R] = polypredcicall_mod(upd,ss2,0.95,'',0.1)
hold on; scatter(upd, ss2,[],'k','filled');
xlabel([rectime ' UPDRS med ON']); ylabel([param ' (s)']);
title([rectime ' med ON vs ' param])

suptitle([rectime ' ' side ' ' condi])
fnm = fullfile(figdir,[upper(rectime) '_' param '_' side '_' condi]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);
end



%--------------------------------------------------------------------------
function perf_updrs_corr(rectime, side,condi,figdir)

global rootdir

load(fullfile(rootdir,'UPDRS_scores.mat'))

SessList = getdata2analyse(rootdir, 'rectype','BEHAV','rectime',rectime, ...
    'patients','allpatients', 'side',side, 'condition',condi);


[~,~,~,~,~,~,perf,stopperf] = calc_RT(SessList,['TrialEvents_nosync_' condi '.mat'],false);

pats = {SessList.patient};
tabpats = updrs_tab.Properties.RowNames;
[~,ai,bi] = intersect(pats,tabpats);

for pk = 1:2
    switch pk
        case 1; % perf
            prf = cell2mat(perf(ai));
            tit = 'perf'
        case 2; %stop perf
            prf = cell2mat(stopperf(ai));
            tit = 'stopperf'
    end
    nns = ~isnan(prf);
    
    fig = figure;
    
    %  MED OFF
    subplot(121)
    upd = updrs_tab{bi,[rectime '_medOFF']};
    % upd = updrs_tab{bi,['preop_medOFF']};
    % [p, R] = polypredcicall_mod(upd(nns),prf(nns),0.95,'robust',0.1);
    [p, R] = polypredcicall_mod(upd(nns),prf(nns),0.95,'',0.1)
    hold on; scatter(upd(nns), prf(nns),[],'k','filled');
    xlabel([rectime ' UPDRS med OFF']); ylabel('performance (%)');
    title([rectime ' med OFF vs ' tit])
    
    
    %  MED ON
    subplot(122)
    upd = updrs_tab{bi,[rectime '_medON']};
    % upd = updrs_tab{bi,['preop_medON']};
    % [p, R] = polypredcicall_mod(upd(nns),prf(nns),0.95,'robust',0.1);
    [p, R] = polypredcicall_mod(upd(nns),prf(nns),0.95,'',0.1)
    hold on; scatter(upd(nns), prf(nns),[],'k','filled');
    xlabel([rectime ' UPDRS med ON']); ylabel('performance (%)');
    title([rectime ' med ON vs ' tit])
    suptitle([rectime ' ' side ' ' condi])
    
    fnm = fullfile(figdir,[upper(rectime) '_' tit '_' side '_' condi '_preop_UPDRS']);
    saveas(fig,[fnm '.jpg'])
    saveas(fig,[fnm '.fig'])
    saveas(fig,[fnm '.pdf'])
    close(fig);
end
end

