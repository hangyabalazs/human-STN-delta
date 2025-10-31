function wcoh_fig(wc,phL,t,f,coi,varargin)
%WCOH_FIG Wavelet coherence map.
%   WCOH_FIG(wc,phL,t,f,coi) draws wav-coherence map.
% 
% Input parameters:
%   WC      MxN matrix of MSWC values. 1D: frequency, 2D: time
%   phL     MxN matrix of MSWC phase lag values. 1D: frequency, 2D: time
%   t       1xN time vector
%   f       1xM vector of freq. components
%   coi     1xN vector for cone of influence

% Johanna Petra Szabó, 05.2025
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

if isempty(varargin);
    cLim = [-.4 .4];
else
    cLim = varargin{1};
end

% fig = figure;
% set(fig,'Visible','on')
h = pcolor(t,log2(f),wc);
colormap(pink)
h.EdgeColor = 'none';
ax = gca;
ytick=round(pow2(ax.YTick),3);
ax.YTickLabel=ytick;
xlabel('Time(s)')
ylabel('Freq (Hz)')

hcol = colorbar;
hcol.Label.String = 'Norm. MSWC';
hold on;
if ~isempty(coi)
    plot(t,log2(coi),'w--','linewidth',2)
end

% if ~isempty(phL)
%     plot_wcoherence_phaselag(phL,wc,3,t,f)
% end

% caxis([0.1 0.4])
caxis(cLim)
end