function SSRTime_PD(sess2analyse)
%SSRTIME_PD estimates SSDp0.5 and SSRT
%   SSRTIME_PD(sess2analyse, ssrt_figdir) calculates the stop signal delay (SSD)
%   where the estimated probability of stopping is 50% probability (SSDp0.5)
%   by fitting a generalized linear regression model to the SSDs corresponding 
%   to successful and failed stop trials. Plots SSDs with fitted model. 
%   Calculates Stop signal reactim time (SSRT) values as the difference between 
%   median reactim time and SSDp0.5.
%   SSDp0.5 and SSRT values are saved to [figdir_pd \ 'Behav' \ 'SSRT' \ [rectime '_' tag]]
%   as 'SSRT_results.mat' for each sessions specified in SESS2ANALYSE struct (see getdata2analyse).
%   
%   See also: GETDATA2ANALYSE, PD_NOSYNCTE

% Johanna Petra Szabó, Panna Hegedus, Balázs Hangya, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

global figdir_pd
rectime = sess2analyse(1).rectime;




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
    
        
    
    % Estimation of SSRT using GLMFIT
    try
        [ssd05 sl] = SSRTime_glmfit(TE);
    catch
        fprintf('NO SSRT estimate possible %s %s %s\n',patnm, side, tag)
        ssd05 = NaN; ssrt = NaN;
    end
    
    if ~isnan(ssd05)
        fig = gcf;
        
        title(patnm)
        
        % Find RT index
        
        rt_patinx =  find(strcmp(patnm,{All_RTs.patient}));
        rt_inx = find(strcmp(side,{All_RTs(rt_patinx).side}));
        
        
        medRT = nanmedian(rmoutliers(All_RTs(rt_patinx(rt_inx)).([rectime '_' tag])));
        meanRT = nanmean(rmoutliers(All_RTs(rt_patinx(rt_inx)).([rectime '_' tag])));
 
        
        ssrt = medRT*1000-ssd05;
        
        
%         hold on; yL = ylim;
%        Ln = line([ssrt ssrt],yL,'Color','k','LineStyle','--')
        
        xL = xlim; yL = ylim;
        t1 = text(xL(2)*.7,yL(2)*.9,['median RT = ' num2str(medRT*1000) ' ms']);
        t2 = text(xL(2)*.7,yL(2)*.8,['mean RT = ' num2str(meanRT*1000) ' ms'] );
        t3 = text(xL(2)*.7,yL(2)*.7,['ssdp0.5 = ' num2str(ssd05) ' ms'] );
        t4 = text(xL(2)*.7,yL(2)*.6,['ssrt = ' num2str(ssrt) ' ms']);
        
%         legend(Ln,'SSRT')
        saveas(fig,fullfile(ssrt_figdir,[patnm,'_',sessnm '_ssd05.png']))
        savefig(fig,fullfile(ssrt_figdir,[patnm,'_',sessnm '_ssd05']))
        close(fig)
    end
    
    if strcmp(rectime,'postop')
        SSRT_results.(patnm).(['S' sessnm]).ssd05 = ssd05;
        SSRT_results.(patnm).(['S' sessnm]).ssrt = ssrt;
    else
        SSRT_results.(patnm).(['S' sessnm]).ssd05 = ssd05;
        SSRT_results.(patnm).(['S' sessnm]).ssrt = ssrt;
    end
    ssrt = [];
    
    
    save(fullfile(ssrt_figdir,'SSRT_results.mat'),'SSRT_results')
    
end



%--------------------------------------------------------------------------
function [ssd05 sig] = SSRTime_glmfit(TE)
%SSRTIME_GLMFIT
%  [ssd05 sig] = SSRTime_glmfit(TE) calculates the stop signal delay (SSD)
%   where the estimated probability of stopping is 50% probability (SSDp0.5)
%   by fitting a generalized linear regression model to the SSDs corresponding 
%   to successful and failed stop trials. Plots SSDs with fitted model. 


% Required input:
%       TE       TrialEvents struct containing behavioral results (see PD_nosyncTE)



x = (TE.StopSignal(~isnan(TE.StopSignal))-0.5)*1000;

Success = intersect(find(TE.SuccesfulStopTrial==1),find(~isnan(TE.StopSignal)));
Fail = intersect(find(TE.FailedStopTrial==1),find(~isnan(TE.StopSignal)));

y = nan(1,length(TE.StimulusOn));
y(Success) = 1; y(Fail) = 0; y = y(~isnan(y));


scatter(x,y,[],'k','filled'); hold on;

[x,rmi ]= rmoutliers(x);
y = y(~rmi);

if ~ismember(0,y)
    yfit = y;
    ssd05 = max(x);
    sig = 0;
else
    
    
    
    b = glmfit(x',y','binomial','link','logit');
    yfit = glmval(b,x','logit');
    
    [sortx sortxi] = sort(x);
    yfit_sort = yfit(sortxi);
    plot(sortx,yfit_sort,'k'); hold on;
    
    
    if yfit_sort(1)>yfit_sort(end)
        sig = 1;
    else
        sig = -1;
    end
    
    
    
    %     try
    [uy uyi] = unique(yfit);
    if max(uy)<0.5 | min(uy)>0.5
        ssd05 = interp1(uy,x(uyi),0.5,'linear','extrap');
    else
        ssd05 = interp1(uy,x(uyi),0.5);
    end
    
    %     catch
    %         ssd05 = (max(x(y==1))+min(x(y==0)))/2;
    %     end
    
    %     if isnan(ssd05) && isempty(find(yfit<0.5))
    %         [miny minyi] = min(yfit);
    %         ssd05 = x(minyi);
    %     end
    
end

scatter(ssd05,0.5,35,'r','filled');
xlabel('Stop Signal Delay (ms)');
ylabel('Probability of stopping');

