function cellnr_pie(X,colors, labels)
%CELLNR_PIE
% CELLNR_PIE(X,colors, labels) creates pie chart with unit nrs stored in X.

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


P = pie(X);
patchHand = findobj(P, 'Type', 'Patch');

if length(X)>2 && size(colors,1)<length(X)
    newColors = [colors; [0.6, 0.6, 0.6]];
    leglabels = {['n= ' num2str(X(1)) ' - ' labels{1}],...
        ['n= ' num2str(X(2)) '-' labels{2}],...
        ['n= ' num2str(X(3)) '-' labels{3}] };
else
    newColors = colors;
    leglabels = arrayfun(@(y) ['n= ' num2str(X(y)) ', ' labels{y}],1:length(X),'UniformOutput',0)
    
end

if ~isempty(colors)
set(patchHand, {'FaceColor'}, mat2cell(newColors, ones(size(newColors,1),1), 3))
end


lgd = legend(leglabels,'Location','bestoutside');