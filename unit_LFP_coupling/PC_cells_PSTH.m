function PC_cells_PSTH(event, fr_name, chanmean, rectype, dominantfreq,PC_win, ds, plot_win)
%PC_CELLS_PSTH  PSTH of sign. phase-coupled units.
% PC_CELLS_PSTH(...) draws 2 types of figures: 
%   1. Maps of the PSTHs of all sign  phase-coupled units, sorted according to 
%       -the magnitude of firing rate change in resp. to the event (ascending order)
%       - coupling strength (MRL )
%       - preferred coupling phase (FTM)
%   2. Average PSTH  
%
% See also: GET_PHAS, SPIK_PHAS_EXTRACTION, PC_CELL_LEVEL
%
% Johanna Petra Szabó, 05.2025
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

[PC_results_dir, mrl, phas, cellids,frnm,varargout] = get_phas(event, fr_name, chanmean, rectype,'signPC', dominantfreq,PC_win, ds, plot_win);
figdir = fullfile(PC_results_dir,'signPC_PSTH');  
if ~isdir(figdir); mkdir(figdir); end;

switch event
    case 'StimulusOn'
        colors = [0 0 1; 0.3010, 0.7450, 0.9330];
        bwin = [-2.5 -1];
    case 'StopSignal'
        colors = [0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250];
        bwin = [-3 -1.5];
    case 'KeyPress1'
        colors = [0.4940, 0.1840, 0.5560; 0.75, 0, 0.75];
        bwin = [-3 -1.5];
    case 'Feedback'
        colors = [0 1 0; 0.4660, 0.6740, 0.1880 ];
        bwin = [1.5 3];
end

for j = 1:3
    switch j;
        case 1; sortix = []; tit2 = 'sorted according to firing rate change'; lab = 'FRsort';
        case 2; sortix = mrl{1}; tit2 = 'sorted according to coupling strength (MRL)'; lab = 'MRLsort';
        case 3; sortix = phas{1}; sortix = mod(sortix,2*pi); tit2 = 'sorted according to pref. phase angle (FTM)'; lab = 'FTMsort';
    end;
    
    f1 = figure;
    [psth] = norm_psth_map1(cellids{1},'activation',event,...
        'baseline','indiv','basl_psth',{},'bwin',bwin, 'parts','all','cLim',[-6 6],...
        'bindex',sortix,'isfig',true);
    setmyplot_balazs(gca);
    
    
    title({['Sign. coupled units around' event ], tit2});
    
    
    fnm = fullfile(figdir, ['PSTH_map_' lab]);
    saveas(f1, [fnm '.jpg'])
    saveas(f1, [fnm '.fig'])
    close(f1);
    
end


f2 = figure;
avg_psth_plot(psth,[], event, colors)
setmyplot_balazs(gca);

title({['Sign. coupled units around' event ], tit2});


fnm = fullfile(figdir, 'PSTH_avg');
saveas(f2, [fnm '.jpg'])
saveas(f2, [fnm '.fig'])
close(f2);

end