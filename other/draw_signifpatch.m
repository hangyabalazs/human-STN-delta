function draw_signifpatch(time,sign_p,col)
%DRAW_SIGNIFPATCH 
% DRAW_SIGNIFPATCH(time,sign_p,col) draws a filled rectangle with reduced opacity
% where significant differences were detected across time.

% Input parameters:
%   TIME    time vector of plotted figure
%
%   SIGN_P  significance vector of 0s and 1s (if 1, significant) with same length as TIME
%
%   COL     color of rectangle
%
% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

yL = ylim;
df = diff([0; sign_p; 0]);

starts = find(df==1);
ends = find(df==-1);
knr = sum(df==1);
if knr==0
    return;
end
if ends(end)>length(time); ends(end) = length(time); end;

hold on;
for k = 1:knr
    patch(time([starts(k) ends(k) ends(k) starts(k)]),  [yL(1)  yL(1)  yL(end)  yL(end)],...
        col,'FaceAlpha',.6,'EdgeColor','none');
end
hold off;