function set_meanFR
% SET_MEANFR   Mean firing rate
%   SET_MEANFR stores mean firing rate property of all detected units in CellBase

% Balázs Hangya, Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

loadcb
cellids = findcell;

propname = 'mean_FR';
if ~ismember(propname,listtag('prop'))
    insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname);
end

for ic = 1:length(cellids)
    act_cellid = cellids{ic};
    
    spikes_file = cellid2fnames(act_cellid,'Spikes');
[pnm,fnm,ext] = fileparts(spikes_file);
if contains(fnm,'TT')&&~contains(fnm,'Ch')
    t = strfind(fnm,'TT');
    fnm(t:t+1) = 'Ch';
    spikes_file = fullfile(pnm,[fnm,ext]);
end

    load(spikes_file)
    TS_sec = TS/10000;
    
    
    secL = TS_sec(end)-TS_sec(1);
    meanFR = length(TS_sec)/secL;
    
    st = setvalue(act_cellid,propname,meanFR); 

end
end