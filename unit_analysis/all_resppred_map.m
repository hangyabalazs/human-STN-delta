function all_resppred_map(EventTypes_all,EventTypes)
% ALL_RESPPRED_MAP    Maps units with all types of behaviorally relevant activity
%   ALL_RESPPRED_MAP(EventTypes_all,EventTypes) draws a map where each type
%   of responsive/ predictive unit is presented (rows). Each column represents a 
%   type of behaviorally relevant activity. Yellow marks activation/ F>S
%   type predictive, blue marks inhibition/ F<S type predictive activity of
%   the unit. Each unit is included only once, but can be assigned with
%   multiple types of activity. 
%
% Input parameters:
%     EVENTTYPES_ALL    1xN cell array of event label for responsive units 
%
%     EVENTTYPES        1xN cell array of event labels for predictive units
%

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu



global cell_dir group_dir

resdir = fullfile(cell_dir,'grouped2','PROPMAPS');
if  ~isfolder(resdir); mkdir(resdir); end;

load(fullfile(group_dir,'RespCells.mat'));
load(fullfile(group_dir,'PredCells.mat'));

prop_titles = {'Resp','Pred'};

propinfo.EventTypes = EventTypes_all; propinfo.caxis = [-1 1];
propinfoS{1} = propinfo;

propinfo.EventTypes = EventTypes; propinfo.caxis = [-1 1];
propinfoS{2} = propinfo;

cells2plot = 'all resp-pred';


pdcellids = cell(4,length(EventTypes_all));
for ei = 1:length(EventTypes_all)
    respevent = EventTypes_all{ei};
    
    
    pdcellids{1,ei} = RespCells.(respevent).none.Activ;
    pdcellids{2,ei} = RespCells.(respevent).none.Inhib;
    if strcmp(respevent,'StopSignal')
        pdcellids{3,ei} = PredCells.([respevent '_StopPartition']).none.F_smaller_than_S;
        pdcellids{4,ei} = PredCells.([respevent '_StopPartition']).none.F_bigger_than_S;
    end
    
end
cellids = unique(cat(2,pdcellids{:}));

draw_propmap3(cellids,cells2plot,prop_titles,propinfoS,resdir);

end



%--------------------------------------------------------------------------
function draw_propmap3(cellids,cells2find,prop_titles,propinfoS,resdir)
propnr = length(prop_titles);
propval = cell(1,propnr);
proplabel = cell(1,propnr);

for k = 1:propnr
    [propval{k}, proplabel{k}] = get_props2plot(prop_titles{k},cellids,propinfoS{k});
end

evsort_inx = contains(proplabel{contains(prop_titles,'Resp')},'KeyPress1'); % sort based on resp to keypress
if ~any(evsort_inx); evsort_inx = 1; end;
if strcmp(cells2find,'all resp-pred') 
    respprop = ismember(lower(prop_titles),'resp'); 
    sortix = [find(propval{respprop}(:,evsort_inx)==1); find(propval{respprop}(:,evsort_inx)==-1); find(propval{respprop}(:,evsort_inx)==0)];
    propval = cellfun(@(x) x(sortix,:),propval,'UniformOutput',0);
end
%%
fig = figure;
for k = 1:propnr
    pnr = length(proplabel{k});
    subplot(3,propnr,[k propnr+k])
    
    imagesc(propval{k})
    caxis(propinfoS{k}.caxis)
    colorbar
        
    if k==1;
        yticks( 1:length(cellids) ); yticklabels(cellids);
    end
    
    xticks(1:pnr); xticklabels( arrayfun( @num2str,1:pnr,'UniformOutput',0 ) );
    
    subplot(3,propnr,[2*propnr+k]); axis off;
    
    %%
    for j = 1:pnr
        text(0,(1-j/pnr),[num2str(j) ' = ' proplabel{k}{j}])
    end
    
    if k==propnr
        subplot(3,propnr,[k propnr+k])
        if any(contains(cat(2,proplabel{:}),'FTM'))
            colormap(hsv);
        end
    end
end
suptitle([cells2find ' cells'] )

set(fig,'Position',get(0,'Screensize'));
savedfnm = dir(fullfile(resdir,['*' cells2find '*']));
if ~isempty(savedfnm)
    try
        fnm = [savedfnm.name(1:end-4) '1'];
    catch
        fnm = [savedfnm(end).name(1:end-4) '1'];
    end
else
    fnm = [cells2find];
end
saveas(fig,fullfile(resdir,[fnm '.jpg']))
saveas(fig,fullfile(resdir,[fnm '.fig']))
close(fig)
end




%--------------------------------------------------------------------------
function [propval, proplabel] = get_props2plot(prop2plot,cellids,propinfo)



% 'resp' | 'pred' | 'PC'
if ~isfield(propinfo,'EventTypes'); propinfo.EventTypes = {'StimulusOn'}; end;



% 'STN_loc'
if ~isfield(propinfo,'IN'); propinfo.IN = true;  end;
if ~isfield(propinfo,'Limit'); propinfo.Limit = true; end;


% PC

if ~isfield(propinfo,'PCprop'); propinfo.PCprop = 'FTM';
elseif ~iscell(propinfo.PCprop); propinfo.PCprop = {propinfo.PCprop}; end;

if ~isfield(propinfo,'win'); propinfo.win = 1;  end;
if ~isfield(propinfo,'downsamp'); propinfo.downsamp = 'no';  end;

if ~isfield(propinfo,'rectype'); propinfo.rectype = {'LFP'};
elseif ~iscell(propinfo.rectype); propinfo.rectype = {propinfo.rectype};end;

if ~isfield(propinfo,'fr_band'); propinfo.fr_band = 'dom_high_delta';  end;
if ~isfield(propinfo,'chanmean'); propinfo.chanmean = true; end;


% Get prop
if contains(prop2plot,'PC')
    pnr = length(propinfo.PCprop);
    rnr = length(propinfo.rectype);
else
    pnr = 1; rnr = 1;
end

enr = length(propinfo.EventTypes);



propval = [];
j = 1;
for rk = 1:rnr
    
    rectype = propinfo.rectype{rk};
    
    for pk = 1:pnr
        if contains(prop2plot,'PC')
            prop2plot2 = [prop2plot propinfo.PCprop{pk}];
        else
            prop2plot2 = prop2plot;
        end
        
        
        for ei = 1:enr
            
            if contains(prop2plot,'PC')
                proplabel{j} = ['PC ' rectype ' ' propinfo.PCprop{pk} ' ' propinfo.EventTypes{ei}];
            elseif contains(lower(prop2plot),'resp') || contains(lower(prop2plot),'pred')
                proplabel{j} = [prop2plot ' ' propinfo.EventTypes{ei}];
            else
                proplabel{1} = prop2plot;
            end
            
            props = get_prop(prop2plot2,cellids,'Events2Align',propinfo.EventTypes(ei),'win',propinfo.win,...
                'downsamp',propinfo.downsamp,'rectype',rectype,'fr_band',propinfo.fr_band,'chanmean',propinfo.chanmean);
            propval = cat(2,propval,props);
            j  = j+1;
        end
        
        
    end
end

end
