function autocorr_PD(cellids,fr_names,freq_bands,win)
%AUTOCORR_PD Auto-correlation
%   AUTOCORR_PD(cellids,fr_names,freq_bands,win) calculates and plots auto-correlations
%       at 0.5 ms resolution for a given window size (WIN) for CELLIDS units.
%       Calculates burst index and oscillation indeces (for details see ACG_MOD).
%       Results are saved in cell_dir\ACG folder.
%       Index values are also saved into CellBase.
%
% Required inputs:
%     CELLIDS       cell array, list of unit IDs
%
%     FR_NAMES      cell array, labels of frequency bands of interest
%
%     FREQ_BANDS    m x 2 matrix, m is the number of freq. bands of interest,
%                   each row corresponds to limits of corresponding frequency bands
%
%     WIN           window size of time period around units to calculate auto-correlation
%
% See also: ACG_MOD, XCORR

% Balázs Hangya, Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


global cell_dir
loadcb


%%
for bi = 1:length(fr_names)
    
    bandname = fr_names{bi};
    band = freq_bands(bi,:);
    
    
    resdir = fullfile(cell_dir,'ACG',bandname);
    if ~isdir(resdir); mkdir(resdir); end;
    
    
    
    
    % acg.m modified to deal with multiple frequency bands
    acg_mod(cellids,win,'issave',true,'resdir',resdir,'band',band,'bandname',bandname,'minspikeno',1000);
    
    
    propname_band = [bandname 'Index'];
    if ~ismember(propname_band,listtag('prop'))
        insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname_band);
    end
    
    load(fullfile(resdir,'ACG_matrices_'));
    for ic = 1:length(cellids)
        st = setvalue(cellids{ic},propname_band,BandIndex(ic));
    end
end


propname_burst = 'BursIndex';
if ~ismember(propname_burst,listtag('prop'))
    insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname_burst);
end

load(fullfile(cell_dir,'ACG','high_delta','ACG_matrices_')); % burst index is the same for each freq band
for ic = 1:length(cellids)
    st = setvalue(cellids{ic},propname_burst,BurstIndex(ic));
end
end
