function component_PC_RT_correlation
% COMPONENT_PC_RT_CORRELATION Finds clusters using spike-phase-coupling strength (MRL) and reaction time value pairs.
%   COMPONENT_PC_RT_CORRELATION fits mix. model of Gaussian distributions
%   on delta spike-phase-coupling strength (MRL) and reaction time value
%   pairs and calculates correlation within clusters. 
%   Plots MRL/ RT value distribution for each cluster. 
%
%
% Johanna Petra Szabó, 10.2025
% Clinic for Neurosurgery and Neurointervention
% Semmelweis University, Budapest, Hungary
% szabo.johanna.petra@semmelweis.hu
%%
compnr = 1; % any nr between 1:5 (nr of fitted components)
cnr = 1; % any nr between 1:5 (nr of cluster)
correlation_within_comp(comp_vals,compnr,cnr,param)

%%
permcomps = [3 1 2]; % switch for plotting, otherwise empty
compare_between_comps(comp_vals,compnr,'rt',permcomps)
compare_between_comps(comp_vals,compnr,'mrl',permcomps)
end

function correlation_within_comp(comp_vals,compnr,cnr,param)
mrl = comp_vals{compnr,cnr}(:,1);
rtval = comp_vals{compnr,cnr}(:,2);

fig = figure;

[p, R] = polypredcicall_mod(mrl,rtval,0.95,'robust',0.01); hold on;
hold on;
scatter(mrl,rtval,[],'k','filled');
xlabel('MRL'); ylabel([param ' (s)']);
title(['GMM for ' num2str(compnr) ' nr. of components. Comp. nr. ' num2str(cnr)])
end

function compare_between_comps(comp_vals,compnr,value,permcomps)
% M = median(rtval);
% mrlup = mrl(rtval>M);
% mrldo = mrl(rtval<=M);
% mup = [mrlup; nan(mL-length(mrlup),1)];
% mdo = [mrldo; nan(mL-length(mrldo),1)];
% boxplot([mup mdo],{'Slow RT', 'Fast RT'});


val = cell(compnr,1);
for j = 1:compnr
    if contains(lower(value),'rt')
        val{j} = comp_vals{compnr,j}(:,2);
    elseif contains(lower(value),'mrl')
        val{j} = comp_vals{compnr,j}(:,1);
    end
end

mL = max(cellfun(@length,val));
val2 = cellfun(@(x) [x; nan(mL-length(x),1)], val, 'UniformOutput',0);

figure
val_bp = cat(2,val2{:});
if ~isempty(permcomps);
    val_bp = val_bp(:,permcomps);
end
boxplot(val_bp);
ylabel(value); xlabel('comp. nr');
set_my_boxplot(gca);
setmyplot_balazs(gca)
[pval_kw,anovatab,stats] = kruskalwallis(val_bp,{},'off');
[c,means,h,gnames] = multcompare(stats,'display','off');

bpax = gca;
for ic = 1:size(c,1)
    comps{ic} = [num2str(c(ic,1)) '-' num2str(c(ic,2))];
    comp_labs{ic} = {gnames{c(ic,1)}, gnames{c(ic,2)}};
    
    pvalc = c(ic,end);
    if pvalc<0.001
        sign = 3;
    elseif  pvalc<0.01
        sign = 2;
    elseif  pvalc<0.05
        sign = 1;
    else
        sign = 0;
    end
    
    
    if sign~=0
        boxplot_astx(bpax,c(ic,1:2),sign)
    end
end

 
% [pval,~,~] = ranksum(rtv{2},rtv{3});
% yL = ylim;
% if pval<0.05; col = 'r'; else; col = 'k'; end;
% text(0.5,yL(2)*.9,['p = ' num2str(pval)],'Color',col);
end