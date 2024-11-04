function PC_pie(event,plot_win,fr_name,downsamp,PC_win,alpha,rectype,chanmean,subregion,PCdir)
%PC_PIE Pie chart of significantly phase-coupled units
%  PC_PIE(event,plot_win,fr_name,downsamp,PC_win,alpha,rectype,chanmean,subregion,PCdir)
%   draws a pie chart respresenting the percentage of significantly
%   phase-coupled units from all detected units, based on spike-phase
%   coupling measures calculated around EVENT, in PC_WIN windows within
%   PLOT_WIN sec time window. 
%
% Input parameters:
%     EVENT             char. array of event label
% 
%     PLOT_WIN          1x2 vector, time window relative to event timestamp in
%                       sec, to use for SPC calculation
% 
%     FR_NAME           char. array, label of frequency band of interest
% 
%     DOWNSAMP          spike or trial nr. downsampled
%                    'no' | 'spike' 
% 
%     PC_WIN            nr of smaller time windows (time window of plot_win will
%                       be divided to PC_WIN nr of smaller windows to use for
%                       SPC calculation) 
%     ALPHA             numeric, alpha level, above this threshold units
%                       are considered sign. coupled)
% 
%     RECTYPE           recording type to use for analyses ('EEG' | 'LFP')
% 
%     CHANMEAN          if true, channel averaged LFP data is used, if false
%                       LFP data is plotted channel-by-channel
%                       (relevant only for intraop LFP data) 
%     SUBREGION         character array or cell array, uses channel data
%                       derived from listed STN subregions
%                       (relevant only for intraop LFP data)
%          'all' | 'Motor' | 'Limbic' | 'Associative' 
%     PCDIR             path to save results and figures
% 
% See also: PC_CELL_LEVEL


% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


if contains(rectype,'LFP')&& chanmean
    chtit = ['chanmean_' subregion 'STN'];
elseif contains(rectype,'LFP')&& ~chanmean
    chtit = 'by-channel';
elseif strcmp(rectype,'EEG')
    chtit = 'F4';
end



resdir = fullfile(PCdir,chtit,event,event,[fr_name '_' num2str(plot_win)]);
load(fullfile(resdir,['PC_results_ds' downsamp '_' num2str(PC_win) 'win.mat']));

if PC_win==1
    ray_sign = structfun(@(x) x.Ray_P<=alpha,PC_results.Hilb_PC);
elseif C_win==2
    ray_sign1 = structfun(@(x) x.Ray_P(1)<=alpha,PC_results.Hilb_PC);
    ray_sign2 = structfun(@(x) x.Ray_P(2)<=alpha,PC_results.Hilb_PC);
    ray_sign = ray_sign1|ray_sign2;
end
signcell = sum(ray_sign);
nosigncell = sum(~ray_sign);
fig = figure;
pie([signcell nosigncell]);
hold on;
legend({['Sign PC, n = ' num2str(signcell)], ['Not sign. PC, n = ' num2str(nosigncell)]});

fnm = fullfile(resdir,'signPC_pie');
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);
