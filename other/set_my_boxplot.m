function set_my_boxplot(ax)
%SET_MY_BOXPLOT Sets the properties of boxplot figures

lines = findobj(ax, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', [0 0 0]);

lins = findobj(ax,'LineStyle','--');
set(lins,'LineStyle','-');

outL = findobj(ax, 'type', 'line', 'Tag', 'Outliers');
set(outL,'Marker','o');
set(outL,'MarkerEdgeColor',[0 0 0]);
set(outL,'MarkerFaceColor',[0 0 0]);
set(outL,'MarkerSize',4);


LWlines = findobj(ax, 'type', 'line', 'Tag', 'Lower Adjacent Value');
set(LWlines,'LineStyle','none')

UWlines = findobj(ax, 'type', 'line', 'Tag', 'Upper Adjacent Value');
set(UWlines,'LineStyle','none')

h = findobj(ax,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),[0.5 0.5 0.5],'FaceAlpha',.3,'LineStyle','none');
   h(j).LineStyle = 'none';
   h(j). Color  = [1 1 1];
end