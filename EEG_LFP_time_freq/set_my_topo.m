function set_my_topo(fig)
%SET_MY_TOPO
%   SET_MY_TOPO(fig) sets properties of topoplot (marker size, line width, contours)

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu



% Electrode location points
Lns = findobj(fig, 'type', 'line');
Lns(1).MarkerSize = .5;

% Head, ear, nose lines
Lns2 = findobj(fig, 'type', 'line','LineStyle','-');
set(Lns2,'LineWidth',0.05);

% Contour lines
Cs = findobj(fig, 'type', 'contour');
set(Cs,'Visible','off')

