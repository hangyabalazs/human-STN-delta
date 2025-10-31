function L_ratio_ID_distrib(L_thr,ID_thr)
% L_RATIO_ID_DISTRIB    Single- / multi unit sorter
%   L_RATIO_ID_DISTRIB(L_thr,ID_thr) extracts the L ratio and ID distribution values
%   (calculated by MCLust during the unit clustering process) of
%   all detected units. Results are saved as 'L_ratio.mat' and 'ID.mat' in
%   cell_dir \ 'L_ratio_ID'   directory. Units are labeled as single- or
%   multi unit based on the specified L_THR L ratio threshold and ID_THR ID
%   distribution threshold. Draws edcf plots for L ratio/ ID values of SUAs
%   and MUAs. 
%   Counts nr. of interspike intervals (ISI) below 1 ms, draws histogram
%   for SUAs and MUAs. 
%

% Balázs Hangya, Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


setpref('cellbase','cell_pattern','Ch');

% Save L ratio and ID distribution of units
save_L_ratio_ID_distrib


%% Sort units to single-/ multi units based on L ratio and ID distribution
[sua_cellids, sua_nr, suaL,suaID,mua_cellids, mua_nr,muaL,muaID] = sort_units(L_thr,ID_thr);


%%
figure;
subplot(121); boxplot(suaL); set(gca,'Yscale','log')
setmyplot_balazs(gca);set_my_boxplot(gca); title('SUA')
subplot(122); boxplot(muaL)
setmyplot_balazs(gca);set_my_boxplot(gca); title('MUA'); set(gca,'Yscale','log')
suptitle('L-ratio')

figure;
subplot(121)
[Ff,Xx] = ecdf(suaL); plot(Xx,Ff,'Color','k','Linewidth',2)
setmyplot_balazs(gca);title('SUA'); xlabel('L threshold'); ylabel('Percent');

subplot(122)
[Ff,Xx] = ecdf(muaL); plot(Xx,Ff,'Color','k','Linewidth',2)

setmyplot_balazs(gca);title('MUA'); xlabel('L threshold'); ylabel('Percent');


% suaID(suaID==Inf) = NaN;
% figure;
% subplot(121); boxplot(suaID); 
% setmyplot_balazs(gca);set_my_boxplot(gca); title('SUA')
% subplot(122); boxplot(muaID); title('MUA')
% setmyplot_balazs(gca);set_my_boxplot(gca);
% suptitle('Isolation Distance')
figure;
subplot(121)
[Ff,Xx] = ecdf(suaID); plot(Xx,Ff,'Color','k','Linewidth',2)
setmyplot_balazs(gca);title('SUA'); xlabel('ID ratio'); ylabel('Percent');

subplot(122)
[Ff,Xx] = ecdf(muaID); plot(Xx,Ff,'Color','k','Linewidth',2)

setmyplot_balazs(gca);title('MUA'); xlabel('ID ratio'); ylabel('Percent');





%% Violation of refractory period (nr. of ISI < 1 ms)
ISI_thr = 1; % ISI threshold in ms
edges = [0:5:50];
[countS,ISI] = countISI(sua_cellids,ISI_thr);
[countM,ISI] = countISI(mua_cellids,ISI_thr);

figure
H1 = histogram(countM,edges); H1.FaceColor = [0 0 0]; H1.FaceAlpha = .8;
hold on;
H2 = histogram(countS,edges); H2.FaceColor = [.5 .5 .5]; H2.FaceAlpha = .8;
legend({'MUA','SUA'}); xlabel('Nr. of spikes <1 ms'); ylabel('Count')

% Rate of units with low percentage of ISI violation (<0.5%)
[M_spS,M_countS,mean_ISIv_rateS, good_rateS]  = ISI_violation_rate(sua_cellids,countS);
[M_spM,M_countM,mean_ISIv_rateM, good_rateM]  = ISI_violation_rate(mua_cellids,countM);

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
function [sua_cellids, sua_nr, suaL,suaID,mua_cellids, mua_nr,muaL,muaID] = sort_units(L_thr,ID_thr)

%%
global cell_dir

resdir = fullfile(cell_dir, 'L_ratio_ID');
if ~isdir(resdir); mkdir(resdir); end;

load(fullfile(resdir, 'L_ratio.mat'))
load(fullfile(resdir, 'ID.mat'))
sua_inx = L_ratio.Lr<L_thr&ID.ID>ID_thr;
mua_inx = ~sua_inx;
sua_cellids = L_ratio.cellids(sua_inx);
mua_cellids = L_ratio.cellids(~sua_inx);

sua_nr = sum(sua_inx);
mua_nr = sum(mua_inx);

unit_types.SUA = sua_cellids;
unit_types.MUA = mua_cellids;


save(fullfile(resdir, 'unit_types.mat'),'unit_types')

suaL = L_ratio.Lr(sua_inx);
suaID = ID.ID(sua_inx);

muaL = L_ratio.Lr(mua_inx);
muaID = ID.ID(mua_inx);
end



%--------------------------------------------------------------------
function [count,ISI] = countISI(pdcells,ISI_thr)

ISI = cell(length(pdcells),1);
count = nan(length(pdcells),1);
for ci = 1:length(pdcells)
    cellid = pdcells{ci};
    
    % Load spike times
    spk = loadcb(cellid,'Spikes');
    spk = spk*1000;
    
    %% Plot histogram of ln(ISI)
    ISI{ci} = diff(spk);
    count(ci) = sum(ISI{ci}<ISI_thr);
end
end

%--------------------------------------------------------------------
function [M_SL,M_count,M_ISIv_rate, good_rate] = ISI_violation_rate(cellids,count)

SLL = nan(length(cellids),1);
for k = 1:length(cellids)
    act_cellid = cellids{k};
    [~,~,chan,unit] = cellid2tags(act_cellid);
    
    
    % Load  spike data
    SP = loadcb(act_cellid,'Spikes');
    SLL(k) = length(SP);
end


M_SL = mean(SLL);
M_count = mean(count);
M_ISIv_rate =(M_count*100)/M_SL;
perc_perU = (count*100)./SLL;
% perc_perUM = mean(perc_perU);
good_units = sum(perc_perU<=0.5);
good_rate = (good_units*100)/length(cellids);
end