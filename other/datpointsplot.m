function sc = datpointsplot(dat,gooddat,lab,color)
% DATPOINTSPLOT     Dots overlayed on boxplots
%   DATPOINTSPLOT(dat,gooddat,lab,color) adds randomly scattered dots
%   corrsponding to each datapoint. Should be plotted over boxplots.
%
%   Input parameters:
%       DAT         MxN matrix, M: nr of datapoints, N: nr of groups (boxplots)
%       GOODDAT     index of groups to plot
%       LAB         cell array of group labels
%       COLOR       color of datapoint dots


% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


if isempty(gooddat); gooddat = 1:size(dat,2); end;
if isempty(color); color = [0 0 0]; end;

dat2 = dat(:,gooddat);

gnr = size(dat2,2);
datlen = size(dat2,1);
for goi = 1:gnr;
    
    sc(goi) = scatter(ones(1,datlen)*goi+randi([-10 10],1,datlen)*0.01,dat2(:,goi),[],color,'filled'); hold on;
    
end
if ~isempty(lab)
xticks(1:gnr); xticklabels(lab(gooddat));
end