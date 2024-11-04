function [tf,pow2,im] = spectr_fig(pow,f,sr,upF,downF,epoch_win, baseline_win, sec2cut,cLims)
%SPECTR_FIG     Spectrogram basic plot
%   [tf,pow2,im] = spectr_fig(pow,f,sr,upF,downF,epoch_win, baseline_win, sec2cut,cLims)
%       Plots a time-frequency spectrogram.
% 
% Required inputs:
%     POW           time-frequency matrix  (1D: frequency, 2D: time)
% 
%     F             frequency vector
% 
%     SR            sampling rate
% 
%     UPF           upper limit of plotted frequencies
% 
%     DOWNF         lower limit of plotted frequencies
% 
%     EPOCH_WIN     1x2 vector, time window (boundaries of pow matrix), relative to events in sec (ex: [-2 2])
% 
%     BASELINE_WIN 	1. 1x2 vector baseline time window relative to events in sec (ex: [-2 -1])
% 
%                   2. 1xn, n>2 vector of baseline time indeces
% 
%     SEC2CUT       seconds to cut from the plot (should be used if wavelet transform was applied on single trial)
% 
% 
% Optional input:
%     CLIMS         color axis limits; if empty 1-99 percentiles of time-freq map
%                   is used (default: [])
% 
% Output: 
%     TF            plotted time-freq map (1D: frequency, 2D: time)
% 
%     POW2          time-freq map prior to baseline normalization
% 
%     FIG           figure handle
%
%

% Johanna Petra Szabó, Hangya Balázs, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


narginchk(8,9)
if nargin<9
    cLims = [];
end
fRange = size(pow,1);
[~,uppLim] = min(abs(f-upF));
[~,downLim] = min(abs(f-downF));


f2_ind = uppLim:downLim;
fValues = f(uppLim:downLim);
% fValues = round(fValues, 1); % Frequency scale


fPos = 1;
if upF>=50 && downF<=30
    [~,uppGamma] = min(abs(fValues-50));
    fPos = [fPos uppGamma];
end
if upF>=30 && downF<=13
    [~,uppBeta] = min(abs(fValues-30));
    fPos = [fPos uppBeta];
end
if upF>=13 && downF<=8
    [~,uppAlpha] = min(abs(fValues-13));
    fPos = [fPos uppAlpha];
end

if upF>=7 && downF<=4
    [~,uppTheta] = min(abs(fValues-7));
    fPos = [fPos uppTheta];
end

%     [~,uppLim] = min(abs(f-foi(2,i)));
%     fPos = [fPos uppLim];

if upF>=4  && downF<4
    [~,uppDelta] = min(abs(fValues-4));
    [~,midDelta] = min(abs(fValues-2));
    [~,midDelta2] = min(abs(fValues-1));
    fPos = [fPos midDelta2 midDelta uppDelta];
end

if length(unique(fPos))==1
    fPos = [uppLim, downLim];
end

fPos = sort(unique(fPos));
if sec2cut~=0
pow2 = pow(uppLim:downLim,sec2cut*sr+1:end-sec2cut*sr);
else
    pow2 = pow(uppLim:downLim,:);
end
tf = pow2;



% fig = figure; 
% fig = [];
% set(fig, 'Visible', 'off');
    
    winlength_sec = abs(epoch_win(1,1))+ abs(epoch_win(1,2));
    winlength = (winlength_sec-sec2cut*2)*sr;
    evpos = (abs(epoch_win(1,1))-sec2cut)*sr;
    
    si1 = sign(epoch_win(1,1));
    win1 = num2str(si1*abs(epoch_win(1,1))-sec2cut);
    si2 = sign(epoch_win(1,2));
    win2 = num2str(si2*abs(epoch_win(1,2))-sec2cut);
    
    % Baseline correction
    
    if numel(baseline_win)==2
        baseline_win_s = baseline_win*sr;
        bas_win_rel = evpos + baseline_win_s+1;
        bas = median(pow2(:,bas_win_rel(1):bas_win_rel(2)),2); % median of baseline period - for each freq comp
        repbas = repmat(bas,[1 size(pow2,2)]);
        tf = ((pow2-repbas)./repbas)*100;
%         pow_bascorr = pow2./repbas;

        xti = [1 bas_win_rel(1,2) evpos winlength];
        xtila =  {strcat(win1,'sec'), strcat(num2str(baseline_win(1,2)),'sec'), 'Event', strcat(win2,'sec')};
    elseif numel(baseline_win)>2
        bas  = baseline_win(f2_ind);
        repbas = repmat(bas,[1 size(pow2,2)]);
        tf = ((pow2-repbas)./repbas)*100;
        
        xti = [1 evpos winlength];
        xtila =  {strcat(win1,'sec'), 'Event', strcat(win2,'sec')};
        
    
    elseif isempty(baseline_win)
        bas_win_rel = [];
        if evpos>=winlength
            xti = [1 winlength];
            xtila =  {strcat(win1,'sec'), strcat(win2,'sec')};
        else
            xti = [1 evpos winlength];
        xtila =  {strcat(win1,'sec'), 'Event', strcat(win2,'sec')};
        end
    end
    
    % spectrogram
    
    % TF - time-frequency map
    
%     cLims = [-150 150];
    im = imagesc(tf), colormap('jet')
    if isempty(cLims)
        cLims = prctile(tf(:),[1,99]);
    end
    set(gca, 'clim', cLims);
    b_rescaleaxis('Y',fValues)
    setappdata(gca,'scaley',fValues);
    b_zoomset_for_wavelet;
    
    colorbar
    
    set(gca,'XTick',xti)
    set(gca, 'XTickLabel',xtila);
    
    tick = fPos;
    ticklab = fValues(fPos);
    set(gca,'YTick', tick)
    set(gca, 'YTickLabel', ticklab)
    % set(gca, 'YTickLabel', '')
    
   ylabel 'Freqeuncy';