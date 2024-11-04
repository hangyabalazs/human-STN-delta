function [RT,Go_RT,FAlarm_RT,Hit,NoStopTrials,StopTrials,perf,stopperf] = calc_RT(sess2analyse,TEfnm,isfig)
%CALC_RT calculates behavioral parameters. 
% CALC_RT(sess2analyse,TEfnm,isfig) calculates behavioral parameters for each
%       trial of patients specified in SESS2ANALYSE struct (see getdata2analyse.m),
%       based on TrialEvents struct (saved as TEFNM mat file in the sess2analyse.folder).
%       If ISFIG is true, trial-by-trial change in behav. parameters is
%       plotted and saved to [figdir_pd/Behav/RT_SSD_trial_change] folder.

%       Output parameters
%           RT      reaction time
%           Go_RT   reaction time of no-stop trials
%           FAlarm_RT   reaction time of failed stop trials (false alarm)
%           Hit     index of correct no-stop trials
%           NoStopTrials    index of no-stop trials
%           StopTrials      index of stop trials
%           perf        performance (ratio of correct responses from all no-stop trials)
%           stopperf        stop performance (ratio of successful stops from all stop trials)

% See also RT_PERF_COMPARE, GETDATA2ANALYSE

% Johanna Petra Szabó, Panna Heged?s, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

dbstop if error

sessnr = size(sess2analyse,1);
[RT,Go_RT,FAlarm_RT,Hit,NoStopTrials,StopTrials,perf,stopperf,successful_stop] = deal(cell(1,sessnr));


for i = 1:sessnr
    
    patnm = sess2analyse(i).patient;
    sessdir = sess2analyse(i).sessfolder;
    tag = sess2analyse(i).tag;
    rectime = sess2analyse(i).rectime;
    side = sess2analyse(i).side;
    
    try
        if contains(TEfnm,'nosync')
            TE = load(fullfile(sessdir,TEfnm));
        else
            TE = load(fullfile(sessdir,TEfnm));
        end
        if strcmp(fieldnames(TE),{'TE'})
            TE = TE.TE;
        end
    catch
        fprintf('No TE % s %s %s\n',sess2analyse(i).patient,...
            sess2analyse(i).side,sess2analyse(i).tag);
        
        [RT{i},Go_RT{i},FAlarm_RT{i},Hit{i},NoStopTrials{i},StopTrials{i},perf{i},stopperf{i},successful_stop{i}] = deal(NaN);
        continue
    end
    
    
    % Calculate RTs
    RT{i} = TE.KeyPress1 - TE.StimulusOn;
    Hit{i} = ~isnan(TE.Correct);
    NoStopTrials{i} = isnan(TE.RealStopTrial); % includes premature responses
    StopTrials{i} = ~isnan(TE.RealStopTrial);
    
    %     perf{i} = sum(Hit{i})/sum(NoStopTrials{i});
%     perf{i} = (sum(Hit{i}(3:end-2))/sum(NoStopTrials{i}(3:end-2)))*100;
    [~,rmx] = rmoutliers(RT{i});
    
    perf{i} = (sum(Hit{i}(~rmx))/sum(NoStopTrials{i}(~rmx)))*100;
    
  
    successful_stop{i} = ~isnan(TE.SuccesfulStopTrial);
    stopperf{i} = (sum(successful_stop{i})/sum(StopTrials{i}))*100;    
 
    
    
    
    FAlarm_RT{i} = RT{i}(StopTrials{i}); % includes stop trials when KeyPress followed Stop Signal
    
    Go_RT{i} = RT{i}(NoStopTrials{i}); % includes premature responses
    
    RT{i}(1) = NaN;
    RT{i}(end) = NaN;% includes also False Alarm!!
    
    TE.RT = RT{i};
    TE.Go_RT = Go_RT{i};
    TE.FAlarm_RT = FAlarm_RT{i};
    
    %     figure; plot(RT{i});
    %     close(gcf);
    %     save(fullfile(sessdir,TEfnm),'TE');
    
    if isfig
        global figdir_pd
        fig = figure;
       ssd_plot(TE);
        if ~isempty(fig)
            fnm = fullfile(figdir_pd,'Behav','RT_SSD_trial_change',[rectime '_' tag]); if ~isfolder(fnm); mkdir(fnm); end;
            saveas(fig,fullfile(fnm,[patnm '_' side '.jpg']))
            saveas(fig,fullfile(fnm,[patnm '_' side '.fig'])); close(fig)
        end
    end
    
end
end

function fig = ssd_plot(TE)
%SSD_PLOT plots trial-by-trial change of behavioral parameters
%   ssd_plot(TE) represents trial-by-trial change of reaction time and stop signal
%       delay by line plots, correct/ incorrect responses and successful/
%       unsuccessful stops are marked by scatter plots (FIG is the handle of the figure). 
%       Behavioral parameters derive from TrialEvents (TE) struct (see PD_nosyncTE.m, RT_perf_compare.m).


    trnr = length(TE.StimulusOn);
    
    ssdx = find(~isnan(TE.StopSignal));
    ssdy = (TE.StopSignal(ssdx)-0.5)*1000;
    
    if isempty(ssdx)
        fprintf('No stop trials\n'); fig = [];
        return
    end
    
    
    rtx = find(~isnan(TE.RT));
    rty = TE.RT(rtx)*1000;
    
    rtinp = interp1(rtx,rty,1:trnr);
    
    
    corrGox = intersect(find(~isnan(TE.Correct)),find(isnan(TE.StopSignal)));
    corrGoy = rtinp(corrGox);
    
    errGox = intersect(find(~isnan(TE.Error)),find(isnan(TE.StopSignal)));
    errGoy = rtinp(errGox);
    

    
    ssdinp = interp1(ssdx,ssdy,1:trnr);
    
    sstopx =  find(~isnan(TE.SuccesfulStopTrial));
    sstopy = ssdinp(sstopx);
    
    fstopx =  find(~isnan(TE.FailedStopTrial));
    fstopy = ssdinp(fstopx);

    fig = figure;
    plot(ssdx,ssdy,'Color','blue'); hold on;
    plot(rtx,rty,'Color','black','LineStyle','--'); hold on;
    scatter(sstopx,sstopy,25,'green','filled'); hold on;
    scatter(fstopx,fstopy,25,'red','filled'); hold on;
%     scatter(corrGox,corrGoy,25,'cyan',"diamond"); hold on;
    scatter(errGox,errGoy,25,'blue',"diamond",'filled');
    
    xlabel('Trial number')
    ylabel('Time (ms)')
    legend({'SSD','RT','SuccesfulStop','FailedStop'})
    ylim([0.2 5]*1000);
end