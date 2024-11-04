function [cellids1 cellids2 bindex1 bindex2 cond_cellidx] = find_cellidx(condition,cell_dir,cellidx1,cellidx2,pdcells)
%FIND_CELLIDX   Cells matching selection criteria
% [cellids1 cellids2 bindex1 bindex2 cond_cellidx] = find_cellidx(condition,cell_dir,cellidx1,cellidx2,pdcells)
%   Finds cellids matching selection criteria (CONDITION) among responsive
%   units.

% Input parameters:
%   CONDITION	character array of selection criteria
%       'bursting' | 'nonbursting' | 'rhythmic' | 'nonrhythmic' | 'none'
%   
%   CELL_DIR	path to directory where information required for selection
%               criteria is stored 
%   
%   CELLIDX1	indeces of activated cells from all cells
%   
%   CELLIDX2	indeces of inhibited cell from all cells
%   
%   PDCELLS     all cells (non-resp + resp)
% 
% Ouput parameters:
%   CELLIDS1      activated cells matching condition
%   
%   CELLIDS2      inhibited cells matching condition
%   
%   BINDEX1       property values associated with cellids1 (ex bursting index)
%   
%   BINDEX2       property values associated with cellids2 (ex bursting index)
%   
%   COND_CELLIDX	all cells matching condition (non-resp + resp)
%
% See also: AUTOCORR_PD, RESPONSESORTER_PD

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu




numCells = length(pdcells);


if strcmp(condition, 'bursting')
    
    load(fullfile(cell_dir,'ACG','high_delta','ACG_matrices_.mat'));
    cond_cellidx = find(BurstIndex>0.1);
    OrigIndex = BurstIndex;
    
elseif strcmp(condition, 'nonbursting')
    
    load(fullfile(cell_dir,'ACG','high_delta','ACG_matrices_.mat'));
    cond_cellidx = find(BurstIndex<=0.2);
    OrigIndex = BurstIndex;
    
    
    
elseif contains(condition, 'rhythmic') && ~contains(condition,'non')
    
    
    try
        bnm = condition(1:strfind(condition,'rhythmic')-1);
        load(fullfile(cell_dir,'ACG',bnm,'ACG_matrices_.mat'));
    catch
        bnm = [upper(condition(1)) condition(2:strfind(condition,'rhythmic')-1)] ;
        load(fullfile(cell_dir,'ACG',bnm,'ACG_matrices_.mat'));
    end
    
    cond_cellidx = find(abs(BandIndex)>0.1);
    OrigIndex = BandIndex;
    
elseif contains(condition, 'rhythmic') && contains(condition, 'non')
    
    if strcmp(condition,'nonrhythmic')
        rhydirs = dir(fullfile(cell_dir,'ACG'));
        
        cond_cellidx0 = [];
        for ri = 3:length(rhydirs)
            
            load(fullfile(cell_dir,'ACG',rhydirs(ri).name,'ACG_matrices_.mat'));
            cond_cellidx0 = [cond_cellidx0; find(abs(BandIndex)>0.1)];
            
        end
        
        OrigIndex = [];
    else
        try
            bnm = condition(4:strfind(condition,'rhythmic')-1);
            load(fullfile(cell_dir,'ACG',bnm,'ACG_matrices_.mat'));
        catch
            bnm = [upper(condition(1)) condition(5:strfind(condition,'rhythmic')-1)] ;
            load(fullfile(cell_dir,'ACG',bnm,'ACG_matrices_.mat'));
        end
        
        cond_cellidx0 = find(abs(BandIndex)>0.1);
        OrigIndex = BandIndex;
    end
    
    cond_cellidx = find(~ismember(1:numCells,cond_cellidx0));
    
elseif strcmp(condition, 'none')
    cond_cellidx = 1:numCells;
    OrigIndex = [];
    
else
    error('non-existant condition')
end


if ~isempty(cellidx1) && ~isempty(cond_cellidx)
    newcellidx1 = intersect(cellidx1,cond_cellidx,'stable');
    cellids1 = pdcells(newcellidx1);
else
    cellids1 = {};
    newcellidx1 = [];
end


if ~isempty(cellidx2) && ~isempty(cond_cellidx)
    newcellidx2 = intersect(cellidx2,cond_cellidx,'stable');
    cellids2 = pdcells(newcellidx2);
else
    cellids2 = {};
    newcellidx2 = [];
end


if ~isempty(OrigIndex) && ~isempty(newcellidx1)
    bindex1 = OrigIndex(newcellidx1);
else
    bindex1 = [];
end
if ~isempty(OrigIndex) && ~isempty(newcellidx2)
    bindex2 = OrigIndex(newcellidx2);
else
    bindex2 = [];
end

