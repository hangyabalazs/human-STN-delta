function corrmapfig(coefD,pvalD,normtype,times,timinx,f,patnm,condi,side,event,evty,act_chan,figdir)
% CORRMAPFIG  Correlation map figure
%   CORRMAPFIG(coefD,pvalD,normtype,times,timinx,f,patnm,condi,side,event,evty,act_chan,figdir)
%   Generates a figure of correlation map COEFD, contours significant
%   changes based on PVALD (map of 1s and 0s, significant where 1).
%
% Input parameters:
%   COEFD       n x m matrix, orrelation coefficiens: frequency components x time
%   PVALD       n x m matrix, significance map: frequency components x time
%   NORMTYPE    char. array, normalization type (already applied at COEFD)
%   TIMES       1 x m vector, time points
%   TIMINX      1 x p vector, time indeces to plot
%   F           1 x n vector, frequency components
%   PATNM       char. array, patient code
%   CONDNI      char. array, task condition
%   SIDE        char. array, tested side
%   EVENT       char. array, event label
%   EVTY        char. array, subevent label
%   ACT_CHAN    char. array, channel label
%   FIGDIR      path to result directory to save figure

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

psig = pvalD<0.05;


fig = figure;
%         set(fig,'Visible','off');
imagesc(times(timinx),1:length(f),coefD(:,timinx));
% pcolor(times(timinx),1:length(f),coefD(:,timinx)); shading interp;
set(gca,'YDir','reverse'); 
colormap(bone);
a = colorbar; ylabel(a,'R','Rotation',270);
caxis([-.4 .4]);
yti = round(linspace(1,length(f), 5)) ;
yticks(yti); yticklabels( arrayfun(@num2str,f(yti),'UniformOutput',0) );

hold on;
contour(times(timinx),1:length(f),psig(:,timinx),'Color','r');
xlabel('Time (s)');
ylabel('Freq (Hz)');
title({[patnm ', ' condi ', ' side],[event ', ' evty],act_chan});

% Save fig
figdir2 = fullfile(figdir,[event '_' evty]);
if ~isfolder(figdir2); mkdir(figdir2); end;

fnm = fullfile(figdir2,[patnm '_' side '_' condi '_' act_chan '_' normtype]);
saveas(fig,[fnm '.jpg']);
saveas(fig,[fnm '.fig']);
saveas(fig,[fnm '.pdf']);
close(fig);



end