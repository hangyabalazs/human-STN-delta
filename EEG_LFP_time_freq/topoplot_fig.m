function topoplot_fig(alldata,chanlocs,sr,plot_win,topobin,cLim,topodir,maskdata,suptit,savetit)
%TOPOPLOT_FIG   Plot scalp maps
%  TOPOPLOT_FIG(alldata,chanlocs,sr,plot_win,topobin,cLim,topodir,maskdata,suptit,savetit)
%   draws a series of scalp maps (topoplot.m by Delorme A & Makeig S (2004)), 
%   splitting PLOT_WIN (whold time window to plot) to TOPOBIN (time bin to 
%   include wihtin one topoplot). Figures are saved in TOPODIR.
%
% Input parameters:
%     ALLDATA   3D matrix: freqs x time x channels
% 
%     CHANLOCS  channel location structure (EEGLAB format)
% 
%     SR        sampling rate
% 
%     PLOT_WIN  1 x 2 vector, time window corresponding to ALLDATA in seconds 
% 
%     TOPOBIN   integer, integer, time window represented on one topoplot in ms (power averaged within time bin)
% 
%     CLIM      1 x 2 vector, color axis limits
% 
%     TOPODIR   directory to save the figure
% 
%     MASKDATA  mask resulting from statistical testing (1 wheren stat. sign., 0 otherwise)
% 
%     SUPTIT    char. array, suptitle of the figure
% 
%     SAVETIT   char. array, file name to save the figure
%
% See also: TOPOPLOT, ERSP_PLOT_STAT, BOOSTAT_EEGLAB_J, TIME_FREQ_PATIENTS

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


plot_win_ts = [1 diff(plot_win)*sr];
topobin_ts = ceil((topobin/1000)*sr);
wind_nr = ceil(diff(plot_win_ts)/topobin_ts);

if isempty(cLim)
    cLim = prctile(alldata(:),[1 99]);
end
% Make plots for each frequency band


%% TOPOPLOT

fig = figure;
set(fig,'Position',get(0,'Screensize'));

st = 1;
spnr = 1;
for wi=1:wind_nr
    
    en = st+topobin_ts-1;
    
    sp(spnr) = subplot(1,wind_nr,spnr);
    xL = xlim; yL = ylim;
    
    
    chanvector = squeeze(mean(alldata(:,st:min(en,plot_win_ts(2)),:),[1 2]));
    
    
    
    
    topoplot(chanvector, chanlocs,'maplimits',cLim);
    set_my_topo(gcf)

    
    
    if ~isempty(maskdata)
        
        mask_topo = squeeze(any(maskdata(:,st:min(en,plot_win_ts(2)),:),[1 2]));
        hold on; scatter(sp(spnr).Children(1).XData(mask_topo),sp(spnr).Children(1).YData(mask_topo),...
            30,'white','filled')
    end
    
    title([num2str((plot_win(1)*1000+topobin*(wi-1))) '-' num2str((plot_win(1)*1000+topobin*wi)) 'ms'])
    
    if wi==wind_nr;  cbar; end;
    
    st = st+ topobin_ts;
    spnr = spnr+1;
    
end
colormap(jet)
set(0, 'DefaultFigureRenderer', 'painters');
suptitle(suptit)
fnm = fullfile(topodir,savetit);
saveas(fig,[fnm '.png']);
saveas(fig,[fnm '.fig']);
% saveas(fig,[fnm '.emf']);
%     saveas(fig,[fnm '.pdf']);
close(fig);

end