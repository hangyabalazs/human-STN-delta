function boxplot_astx(bpax,compbox,sign)
%BOXPLOT_ASTX adds asterisks above boxplots
% BOXPLOT_ASTX(bpax,compbox,sign) draws asterisks above bar/boxplots where
% results are statistically significant. 
%   Input arguments:
%       BPAX - bar/boxplot axis
%       COMPBOX - index of compared groups (each boxplot corresponds to a group)
%       SIGN - strength of significance (p<0.05 -> sign = 1; p<0.01 -> sign = 2; p<0.01 -> sign = 3)

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

yt = get(bpax, 'YTick');
axis([xlim    0  ceil(max(yt)*1.2)])
xt = get(bpax, 'XTick');
hold on


plot(xt(compbox), [1 1]*max(yt)*1.1, '-k'); hold on;
text(mean(xt(compbox)),max(yt)*1.15, repmat('*',[1 sign]),'Color','k')
hold off