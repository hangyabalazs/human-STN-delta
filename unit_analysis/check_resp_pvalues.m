function check_resp_pvalues(EventTypes)
% CHECK_RESP_PVALUES Q-Q plot of observed and expected p values obtained by testing unit responsiveness.
%   check_resp_pvalues(EventTypes) gets p values resulting from testing
%   unit responsiveness to events included in cell array of EVENTTYPES,
%   plots observed and expected p value pairs as a scatter plot.
%
% See also: RESPONSESORTER_PD
%
% Johanna Petra Szabó, 10.2025
% Clinic for Neurosurgery and Neurointervention
% Semmelweis University, Budapest, Hungary
% szabo.johanna.petra@semmelweis.hu


global stat_dir
pdcells = findcell;

%%
for iE = 1:length(EventTypes)
    event = EventTypes{iE};
    resdir = fullfile(stat_dir,event);
    
    get_pvalues(pdcells,event,resdir)
    
end


end


function get_pvalues(pdcells,event,resdir)

partition = 'all';


numCells = length(pdcells);
[inh_p,act_p,bp_p] = deal(nan(1,numCells));
for iC = 1:numCells
    cellid = pdcells{iC};   % current cell
    
    cellidt = regexprep(cellid,'\.','_');
    try
        
        fnm = fullfile(resdir,[cellidt '_' event '_' partition(2:end) '*.mat']);
        
        curr_fnm = dir(fnm);
        load(fullfile(resdir,curr_fnm.name));
    catch
        continue
    end
    
    close(gcf);
    inh_p(iC) = stats1{1}.Wpi;
    act_p(iC) = stats1{1}.Wpa;
    stats1 = {};
    
    fnm_B = fullfile(resdir,[cellidt '_' event '_' partition(2:end) '*boxstat*.fig']);
    curr_fnmB = dir(fnm_B);
    openfig(fullfile(resdir,curr_fnmB.name));
    tH =findall(gcf,'Type','Text');
    try
        bp_p(iC) = str2num(tH(~cellfun(@isempty,{tH.String})).String);
    catch
        disp('no')
    end
    close(gcf)
end

inh_fdr_alpha = fdr(inh_p,0.05);
act_fdr_alpha = fdr(act_p,0.05);

sign_inh = find(inh_p<0.05);
sign_act = find(act_p<0.05);
intersect(sign_inh, sign_act)


exp_pvals = (1:numCells)' / (numCells + 1);  % Avoids 0 and 1 extremes

figure;
for k = 1:3
    switch k; case 1; pval = inh_p; tit = [event '-inhib. units'];
        case 2; pval = act_p; tit = [event '-activ. units'];
        case 3; pval = bp_p; tit = [event '- MannWhitney'];
    end;
    
    [sorted_pvals, sortix] = sort(pval);
    subplot(1,3,k);
    plot(exp_pvals, sorted_pvals, 'bo')  % Blue circles
    hold on;
    plot([0 1], [0 1], 'k--')  % Reference line y = x
    xlabel('Expected p-values (Uniform)');
    ylabel('Observed p-values');
    title(tit);
    
end
suptitle('Q-Q Plot of Observed vs. Expected p-values')

figure
for k = 1:3
    switch k; case 1; pval = inh_p; tit = [event '-inhib. units'];
        case 2; pval = act_p; tit = [event '-activ. units'];
        case 3; pval = bp_p; tit = [event '- MannWhitney'];
    end;
    
    subplot(1,3,k)
    histogram(exp_pvals,100);
    hold on;
    histogram(pval,100);
    title(tit);
end



end