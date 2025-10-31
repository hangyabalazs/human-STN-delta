function h = stacked_bar_perc(y,barcolors,textcolor)
%STACKED_BAR_PERC   Stacked bar plots with percentage 
%   STACKED_BAR_PERC(y,barcolors,textcolor) draws 2 groups of stacked bar plots using
%   data in Y. Percentage of stacked bars of the stacked bar groups are
%   displayed, aligned to the middle of each stacked bar. 
% 
% Input parameters:
%   Y           data to plot
%   BARCOLORS   2 x 3 color matrix
%   TEXTCOLOR   1 x 3 color vector or color string for text
%
% See also: BAR
%
% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


h = bar(y,'stacked', 'FaceColor','flat');
for k = 1:size(barcolors,1)
    h(k).CData = barcolors(k,:);
end

% Compute percentage
yp = y./sum(y,2) * 100; 
% Compute bar segment centers
xbarCnt = vertcat(h.XEndPoints);
ybarTop = vertcat(h.YEndPoints);
ybarCnt = ybarTop - (y/2)';
% Create text strings
txt = compose('%.1f%%',yp);
% Add text
th = text(xbarCnt(:), ybarCnt(:), reshape(txt',[numel(txt) 1]), ...
    'HorizontalAlignment', 'center', ....
    'VerticalAlignment', 'middle', ...
    'Color', textcolor,....
    'FontSize', 7);