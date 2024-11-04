function PD_ssrt_behav_MAIN(conditions, all_groups)
% PD_SSRT_BEHAV_MAIN Main function for behavioural analysis of Stop Signal
% Reaction Time task of Parkinsonian patients implanted with deep brain
% stimulator. 
%
% Input parameters:
%       CONDITIONS      Nx2 cell array of task conditions to compare
%               -first column corresponds to recording time, second column to
%               DBS stimulation condition
%               -ex.:{'preop','stimoff';'intraop','stimoff';'postop','stimoff';'postop','stimon'};
%       ALL_GROUPS      cell array of patient groups; each cell contains     
%                       the cell array of subgroups within that patient group

% See also: PD_NOSYNCTE, RT_PERF_COMPARE, SSRTIME_PD, SSDP05_COMPARE, PREOP_POSTOP_RT, UPDRS_BEHAV_CORR

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

global filesdir

condnr = size(conditions,1);
stat = 'nonparam'; % statistics (ranksum if 'nonparam', ttest2 if 'param')

%%% Save TrialEvents matrix NOT synchronized  with EEG

for coci = 1:condnr; % loop over conditions
    rectime = conditions{coci,1};
    condi = conditions{coci,2};
    
    sess2analyse = getdata2analyse(filesdir, 'rectype','BEHAV',...
        'rectime',rectime,'patients', 'allpatients', 'side','bothside', 'condition',condi);
    
    for snr = 1:size(sess2analyse,1)
        
        % Dataset parameters
        patnm = sess2analyse(snr).patient;
        side = sess2analyse(snr).side;
        tag = sess2analyse(snr).tag;
        sessdir = sess2analyse(snr).sessfolder; stf = strfind(sessdir,'\');
        sessnm = sessdir(stf(end)+1:end);
        
        
        
        current_bhfile = dir([sessdir '\' 'pdtask*.mat']); % find behavior files
        
        
        if length(current_bhfile)==1
            load(fullfile(sessdir,current_bhfile.name));
        elseif  length(current_bhfile)==2
            tx= find(contains({current_bhfile.name},tag));
            load(fullfile(sessdir,current_bhfile(tx).name));
        else
            fprintf('NO BEHAV %s, %s, %s\n',patnm,sessnm,side)
            continue
        end
        
        TE = PD_nosyncTE(Results,1);
        save(fullfile(sessdir,['TrialEvents_nosync_' tag '.mat']),'TE');
    end
end



%%% Comparison of reaction time and performance within and across patients

remoutliers = 'indiv'; % 'indiv' = remove outlier trials by individual patients| 'avg' = remove outliers patients| 'both'
avgfig = true; % = true to make figures with all patients
addlines = false; % = true to represent patients with lines interconnecting conditions
addpoints = true; % = true to represent patients as points
rtplots = true; % = true to plot RT figures
perfplots = true; % = true to plot performance figures

for  g = 1:length(all_groups)
    groups_nm = all_groups{g};
    
    if contains(groups_nm,'all');
        indivfig = true; % = true to make figures for each patient
    else
        indivfig = false;
    end
    
    RT_perf_compare(groups_nm,conditions,stat,remoutliers, indivfig,avgfig,addlines,addpoints,rtplots,perfplots)
end


%%% Estimate SSD p0.5

for coci = 1:condnr; % loop over conditions
    rectime = conditions{coci,1};
    condi = conditions{coci,2};
    
    sess2analyse = getdata2analyse(filesdir, 'rectype','BEHAV',...
        'rectime',rectime,'patients', 'allpatients', 'side','bothside', 'condition',condi);
    
    SSRTime_PD(sess2analyse)
end


%%% Comparison of SSDp0.5

for  g = 1:length(all_groups)
    groups_nm = all_groups{g};
    
    SSDp05_compare(groups_nm,'ssd05',stat,conditions)
end


%%% RT change relative to preop condition

for  g = 1:length(all_groups)
    groups_nm = all_groups{g};
    if contains(groups_nm,'all'); continue; end;
    
    preop_postop_RT(groups_nm,conditions)
end

%%% Correlate behavioural parameters with UPDRS scores

updrs_conditions = {'preop','stimoff';'postop','stimon'};

updrs_behav_corr(updrs_conditions)
end