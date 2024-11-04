function L_ratio_ID_distrib(L_thr,ID_thr)
% L_RATIO_ID_DISTRIB    Single- / multi unit sorter
%   L_RATIO_ID_DISTRIB(L_thr,ID_thr) extracts the L ratio and ID distribution values 
%   (calculated by MCLust during the unit clustering process). 
%   all detected units. Results are saved as 'L_ratio.mat' and 'ID.mat' in 
%   cell_dir \ 'L_ratio_ID'   directory. Units are labeled as single- or
%   multi unit based on the specified L_THR L ratio threshold and ID_THR ID
%   distribution threshold. 
%   

% Balázs Hangya, Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


setpref('cellbase','cell_pattern','Ch');

% Save L ratio and ID distribution of units
save_L_ratio_ID_distrib


% Sort units to single-/ multi units based on L ratio and ID distribution
[sua_cellids, sua_nr, mua_cellids, mua_nr] = sort_units(L_thr,ID_thr);

end


%--------------------------------------------------------------------------
function save_L_ratio_ID_distrib

global cell_dir

resdir = fullfile(cell_dir, 'L_ratio_ID');
if ~isdir(resdir); mkdir(resdir); end;

cellids = findcell;
feature_names = {'Amplitude','Energy','Time','WavePC1','WavePC2','WavePC3'};

[DM, LrC] = deal( nan(length(cellids),1 ) );
for k = 1:length(cellids)
    cid = cellids{k};
    try
    [DM(k), LrC(k), ~] = LRatio_pd_ssrt(cid,feature_names,1);
    catch
        disp(k)
    end
    
end

fig = figure;
subplot(211)
histogram(DM,[0:5:200 300:100:1000 2000:1000:5000]); title('Isolation Distance');
subplot(212)
histogram(LrC,[0:1/100:0.3 0.4:1/20:1 2:1:9]); title('L-ratio');
saveas(fig,fullfile(resdir,'LR_ID.jpg'))
saveas(fig,fullfile(resdir,'LR_ID.fig'))
saveas(fig,fullfile(resdir,'LR_ID.pdf'))
close(fig);

L_ratio.cellids = cellids;
L_ratio.Lr = LrC;

ID.cellids = cellids;
ID.ID = DM;

save(fullfile(resdir, 'L_ratio.mat'),'L_ratio')
save(fullfile(resdir, 'ID.mat'),'ID')
end



%--------------------------------------------------------------------------
function [sua_cellids, sua_nr, mua_cellids, mua_nr] = sort_units(L_thr,ID_thr)

%%
global cell_dir

resdir = fullfile(cell_dir, 'L_ratio_ID');
if ~isdir(resdir); mkdir(resdir); end;

load(fullfile(resdir, 'L_ratio.mat'))
load(fullfile(resdir, 'ID.mat'))
sua_inx = L_ratio.Lr<L_thr&ID.ID>ID_thr;
sua_cellids = L_ratio.cellids(sua_inx);
mua_cellids = L_ratio.cellids(~sua_inx);

sua_nr = sum(sua_inx);
mua_nr = sum(~sua_inx);

unit_types.SUA = sua_cellids;
unit_types.MUA = mua_cellids;


save(fullfile(resdir, 'unit_types.mat'),'unit_types')

end

