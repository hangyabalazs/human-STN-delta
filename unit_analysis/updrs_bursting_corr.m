function updrs_bursting_corr(avgBI, figdir)
%UPDRS_BEHAV_CORR calculates correlation between UPDRS scores and bursting index
%   UPDRS_BEHAV_CORR calculates Pearson's correalation between preop and postop (DBS-on)
%   UPDRS scores (tested both with and without medication) and bursting
%   index of units. 
%   Creates scatter plots and overlays regression line (robust regression) with confidence interval.
%
%   Input parameters
%       AVGBI   true | false, if true, median of bursting indeces
%               corresponding to each patient is consdiered
%       FIGDIR  result directory to save plots
%
% See also UPDRS_BEHAV_CORR, ACG_MOD, POLYPREDCICALL_MOD

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


for rt = 1:2
    switch rt; case 1; rectime = 'preop'; case 2; rectime = 'postop'; end;
    for oo = 1:2
        switch oo; case 1; offon = 'off'; case 2; offon = 'on'; end;
        bursting_updrs(rectime,offon,figdir,avgBI)
    end
end

end


%--------------------------------------------------------------------------
function bursting_updrs(rectime,offon,figdir,avgBI)

global rootdir

load(fullfile(rootdir,'UPDRS_scores.mat'))
upd = updrs_tab{:,[rectime '_med' upper(offon)]};

allpats = updrs_tab.Properties.RowNames;


patnr = length(allpats);
medbrst = nan(patnr,1);
allbrst = []; allupd = [];
for pati = 1:patnr
    patnm = allpats{pati};
    patcells = findcell('RAT',patnm);
    brst = get_prop('bursting',patcells);
    %     figure;
    %     boxplot(brst,'plotstyle','compact'); hold on; scatter(ones(1,length(brst)),brst,[],'m','filled'); ylim([-0.8 0.8]);
    %     keyboard;
    %     pause(1);
    %     close(gcf);
    % [~,ix] = sort(brst,'ascend');
    % patcell_sort = patcells(ix)';
    
    if avgBI
        bval = nanmedian(brst);
    else
        bval = brst;
    end
    allbrst = [allbrst; bval];
    allupd = [allupd; repmat(upd(pati),length(bval),1)];
    
end

nn = ~isnan(allbrst);
allbrst2 = allbrst(nn);
allupd2 = allupd(nn);

fig = figure;



% [p, R] = polypredcicall_mod(allbrst2,allupd2,0.95,'robust',0.1);
[p, R] = polypredcicall_mod(allbrst2,allupd2,0.95,'',0.1);
hold on; scatter(allbrst2, allupd2,[],'k','filled');
xlabel('bursting index'); ylabel([rectime ' med ' upper(offon) ' UPDRS'])



suptitle([rectime ' med-' offon ])
fnm = fullfile(figdir,[upper(rectime) '_med' upper(offon) '_medianBI' char(string(avgBI))]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);
end
