function get_SSDs_RTs(conditions)
%GET_SSDS_RTS Stop Signal Delay and Reaction Time values for each patient
%
% Input parameters:
%       CONDITIONS      Nx2 cell array of task conditions to compare
%               -first column corresponds to recording time, second column to
%               DBS stimulation condition
%               -ex.:{'preop','stimoff';'intraop','stimoff';'postop','stimoff';'postop','stimon'};
%
% Johanna Petra Szabó, 10.2025
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

global filesdir

condnr = size(conditions,1);

%%
for coci = 1:condnr; % loop over conditions
    rectime = conditions{coci,1};
    condi = conditions{coci,2};
    
    sess2analyse = getdata2analyse(filesdir, 'rectype','BEHAV',...
        'rectime',rectime,'patients', 'allpatients', 'side','right', 'condition',condi);
    
    [SSD_med, SSD_std, RT_med, RT_std] = get_SSDs_RTs_cond(sess2analyse);
    SSD_med = SSD_med/1000;
    SSD_std = SSD_std/1000;
    
end
end


function [SSD_med, SSD_std, RT_med, RT_std] = get_SSDs_RTs_cond(sess2analyse)

global figdir_pd

rectime = sess2analyse(1).rectime;


[SSD_med, SSD_std, RT_med, RT_std] = deal(nan(size(sess2analyse,1),1));

rt_figdir = fullfile(figdir_pd, 'Behav','RTmedian');
load(fullfile(rt_figdir,['All_RTs_all.mat']));


for snr = 1:size(sess2analyse,1)
    
    patnm = sess2analyse(snr).patient;
    side = sess2analyse(snr).side;
    tag = sess2analyse(snr).tag;
    sessdir = sess2analyse(snr).sessfolder; stf = strfind(sessdir,'\');
    sessnm = sessdir(stf(end)+1:end);
    
    ssrt_figdir = fullfile(figdir_pd,'Behav','SSRT',[rectime '_' tag]); if ~isdir(ssrt_figdir); mkdir(ssrt_figdir); end;
    
    
    if exist(fullfile(ssrt_figdir,'SSRT_results.mat')) == 2
        load(fullfile(ssrt_figdir,'SSRT_results.mat'));
    else; SSRT_results = struct; end;
    
    try
        
        TE = load(fullfile(sessdir,['TrialEvents_nosync_' tag '.mat']));
        if strcmp(fieldnames(TE),{'TE'})
            TE = TE.TE;
        end
        
    catch
        fprintf('No TE file %s,%s,%s\n',patnm,side,tag)
        continue
    end
    
    
    SSDs = (TE.StopSignal(~isnan(TE.StopSignal))-0.5)*1000;
    
    
    rt_patinx =  find(strcmp(patnm,{All_RTs.patient}));
    rt_inx = find(strcmp(side,{All_RTs(rt_patinx).side}));
    RTs = rmoutliers(All_RTs(rt_patinx(rt_inx)).([rectime '_' tag]));
    
    SSD_med(snr) = nanmean(SSDs);
    SSD_std(snr) = nanstd(SSDs);
    
    RT_med(snr) = nanmean(RTs);
    RT_std(snr) = nanstd(RTs);

end
end