function phase_hist_sinewave_plot(event,plot_win,PCdir,binnr)
%PHASE_HIST_SINEWAVE_PLOT Population phase histogram + preferred phase with sine wave
%   PHASE_HIST_SINEWAVE_PLOT(event,plot_win,PCdir,binnr) draws the
%   population histogram of preferred phases of significantly coupled
%   units. The preferred phase of the population is marked. A sine wave is
%   overlaid on the histogram. 
% 
% Input parameters:
%     EVENT             char. array of event label
% 
%     PLOT_WIN          1x2 vector, time window relative to event timestamp in
%                       sec, to use for SPC calculation
%     PCDIR             path to save results and figures
% 
%     BINNR             numeric, nr. of bins in histogram
% See also: PC_GROUPS


% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

sinefreq = 1;
times = 0:1/100:1;
wave = sin(2*pi*sinefreq*times);


load(fullfile(PCdir,'by-channel', event, event,['dom_high_delta_' num2str(plot_win)],'Groups_modelselection','PC_groups_ds_no_win1.mat'));
resdir = fullfile(PCdir,'by-channel', event, event,['dom_high_delta_' num2str(plot_win)]);
phasvals = PC_groups.signPC.phasevals{1};

fig = figure;


edges = linspace(0,2*pi,binnr);
phasvals = mod(phasvals,2*pi);


yyaxis left
h = histogram(phasvals,edges,'Normalization','probability');
h.FaceColor = [.5 .5 .5]; set(gca,'YColor','k'); ylabel('Norm. count');
yyaxis right

plot(linspace(0,2*pi,101),wave,'LineWidth',3,'Color','k')
yticklabels({}); yticks([]); set(gca,'YColor','k');
xlim([0 2*pi]); xlabel('Phase');
xticks([0 pi/2 pi 3*pi/2 2*pi]); xticklabels(arrayfun(@num2str,[0 pi/2 pi 3*pi/2 2*pi],'UniformOutput',0));

fnm = fullfile(resdir,['signPC_phasehist_sinewave_plot_' event]);
saveas(fig,[fnm '.fig']);
saveas(fig,[fnm '.jpg']);
saveas(fig,[fnm '.pdf']);
close(fig);