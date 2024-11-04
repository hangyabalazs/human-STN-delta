function varargout = avg_psth_plot(psth_R1,psth_R2, alignevent, colors,varargin)
%AVG_PSTH_PLOT  Average peri-event time histogram
%   AVG_PSTH_PLOT(psth_R1,psth_R2, alignevent, colors,...) draws averaged
%   PSTHs with standard error for data in PSTH_R1 and PSTH_R2.
% 
%   Required inputs:
%       PSTH_R1, PSTH_R2   n x t matrices, n: number of units, t: number of timepoints; 
%                          if not empty, t of PSTH_R1 and PSTH_R2 has to be the same   
%       ALIGNEVENT         char. array of event label
%       COLORS             2 x 3 matrix, first row: color for data in
%                          PSTH_R1, second row: color for data in PSTH_R2
%       
%   Optional input:
%        WN                1x2 vector, time limits in sec; time vector is
%                          calculated with 1000 Hz sampling frequency
%                          (default value: [-3 3])
%   

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu



% if ~isempty(varargin); yL = varargin{1}; else; yL = [-7 7]; end;
if nargin<5
    wn = [-3 3]; % time window
else
    wn = varargin{1};
end

% Time vector 

dt = 0.001; % time resolution
time = wn(1):dt:wn(2);   % time vector




% Plot average PSTH
psth_R1(psth_R1==Inf) = NaN;
psth_R2(psth_R2==Inf) = NaN;
if size(psth_R1,1)>1
   p1 = errorshade(time,nanmean(psth_R1,1),nanstd(psth_R1,1)/sqrt(size(psth_R1,1)),...
        'LineColor',colors(1,:),'ShadeColor',colors(1,:)); hold on % excitation
elseif ~isempty(psth_R1)
    p1 = plot(time,psth_R1,'Color',colors(1,:)); hold on % excitation
elseif isempty(psth_R1)
    p1 = plot([]);
end
if size(psth_R2,1)>1
    p2 = errorshade(time,nanmean(psth_R2,1),nanstd(psth_R2,1)/sqrt(size(psth_R2,1)),...
        'LineColor',colors(2,:),'ShadeColor',colors(2,:)); hold on % inhibition
elseif ~isempty(psth_R2)
    p2 = plot(time,psth_R2,'Color',colors(2,:)); hold on % excitation
elseif isempty(psth_R2)
    p2 = plot([]);
end
set(gca,'XLim',wn);
% set(gca,'YLim',yL);
yL = ylim;
notnan1 = ~any(isnan(psth_R1),2);
notnan2 = ~any(isnan(psth_R2),2);

if ~isempty(psth_R1)
    txt1 = text(wn(2),yL(2)*0.9,['n = ' num2str(sum(notnan1)) '/'  num2str(size(psth_R1,1)) ],'Color',colors(1,:),'HorizontalAlignment','right');
else
    txt1 = [];
end
if ~isempty(psth_R2)
    txt2 = text(wn(2),yL(2)*0.8,['n = ' num2str(sum(notnan2)) '/' num2str(size(psth_R2,1))],'Color',colors(2,:),'HorizontalAlignment','right')
else
    txt2 = [];
end
xlabel(['Time from ' alignevent '(ms)'])
ylabel('Norm. firing rate')



varargout{1} = p1;
varargout{2} = p2;
varargout{3} = txt1;
varargout{4} = txt2;