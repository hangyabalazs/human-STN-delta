function [props, cellids2] = get_prop(prop,cellids,varargin)
%GET_PROP finds required property values of selected units
%   [props, cellids2] = GET_PROP(prop,cellids,varargin)
% Gets specified properties (prop) of cellids
% In case of phase-coupling, peak amplitude and - latency properties
% cell array of event labels (EventTypes) should also be specified.
%   Required input parameters:
%       PROP    character array, label of unit property
%           - bursting index: 'bursting'  
%           - spike-phase coupling: 'PC MRL' (mean resultant length)| 
%                   'PC FTM' (first trigonometric moment)| 'PC Pval' (p value, Reighley's test)| 
%                   'PC sign' (if significantly coupled)| 
%           - localization within the STN: 'STN_loc' 
%           - responsive/ predictive units: 'resp'  | 'pred' | 'latency' | 'amplitude'
%           - firing rate: 'rate' 
%
%       CELLIDS     n x 1 cell array of unit IDs; 
%                   n x 2 cell array, when PROP is 'latency' | 'amplitude' | 'rate', 
%                   first colomn corresponds to activation/max value,
%                   second column corresponds to inhibition/min value.
%
%
%   Optional input parameters (name-value pairs with default values)
%
%               'Events2Align'      cell array of event(s) used to align data;
%                           	 in case of properties related to
%                               phase-coupling and responsive/ predictive units
%                               (default value: {'StimulusOn'})
%
%        Relevant only when prop. is phase-coupling-related:
%
%               'win'       if 1, get PC calculated in whole data epoch (-2 2 s around event) 
%                           if 2, get PC calculated separately 
%                           in time windows before (-2 0) and after (0 2) event 
%                           (default value: 1)
% 
%               'downsamp'  spike or trial nr. downsampled  
%                    'no' | 'spike' (default value: 'no')
% 
%               'rectype'   recording type
%                    'EEG' | 'LFP' ('default value: 'LFP')
% 
%               'fr_band'   frequency band (default value: 'delta')
%
%               'chanmean'  if true, PC derived from averaged data of LFP channels
%                           if false, PC derived from individual LFP
%                           channels
%                    true | false (default: true)
%
%        Relevant only when prop. is responsive/predictive unit-related:
% 
%               'sigma'     smoothing parameter of PSTHs (default value: 0.01)
%
%               'test_win'   1 x 2 vector, time window around event (sec)
%                           used for testing responsiveness of units
%                           (default value: [0 1])
%
%               'resptype'   response type
%                    'activation' | 'inhibition' (default value: 'activation')
% 
% 
%   Output parameters:
%       PROPS       - m x p matrix, where m is the number of units with available
%                   property value; p is the number of events defined in 'Events2Align'
%                   - 2 x 1 cell array when CELLIDS is n x 2 cell array
%
%       CELLIDS2    corresponding cellids with available property value
% 
%
% See also: spike_phase_coupling_PD, responsesorter_PD


% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu



prs = inputParser;
addRequired(prs,'prop',@ischar)
addRequired(prs,'cellids',@iscell)
addParameter(prs,'Events2Align',{'StimulusOn'},@iscell)
addParameter(prs,'win',2,@isnumeric)
addParameter(prs,'downsamp','no',@ischar)
addParameter(prs,'rectype','LFP',@ischar)
addParameter(prs,'fr_band','delta',@ischar)
addParameter(prs,'chanmean',true,@islogical)
addParameter(prs,'sigma',0.01,@isnumeric)
addParameter(prs,'test_win',[0 1],@isvector);
addParameter(prs,'resptype',{'activation'},@(x) iscell(x));
parse(prs,prop,cellids,varargin{:})
p = prs.Results;


global rootdir cell_dir

cellids2 = p.cellids;

if strcmp(lower(p.prop),'bursting')
    props = getvalue('BursIndex',p.cellids);
    
    
elseif contains(p.prop,'PC');
    
    props = [];
    for evi = 1:size(p.Events2Align,2)
        [MRLs, Ps, FTM] = get_PCvalues(p.cellids,p.Events2Align(:,evi),p.fr_band,p.rectype,p.win,p.downsamp,p.chanmean);
        
        if contains(p.prop,'MRL')
            px = permute( MRLs,[2 1] );
        elseif contains(p.prop,'Pval')
            px = permute( Ps,[2 1] );
            
        elseif contains(p.prop,'FTM')
            px = permute( FTM,[2 1] );
        end
        
        if contains(p.prop,'sign')
            px(Ps>0.05) = NaN;
        end
        
        props = cat(2,props,px);
    end
    
    
elseif contains(p.prop,'STN_loc')
    
     load(fullfile(rootdir,'STN_loc','STN_subregions.mat'));
     load(fullfile(rootdir,'STN_loc','STN_loc.mat'));
     subreg_names = STN_loc.subreg_names;
     props = nan(length(cellids),1);
     for ic = 1:length(cellids)
         cid = cellids{ic};
         [patnm,sessnm,chan,~] = cellid2tags(cid);
         if strcmp(sessnm(end),'l'); side = 'left'; else; side = 'right'; end;
         dist = STN_subregions.Patients.(patnm).(side).(['Ch' num2str(chan)]).Dist;
         
         for k = 1:length(subreg_names)
             csin = STN_loc.SubReg_Cellids.(subreg_names{k}).IN; % in subregion
             cslim = STN_loc.SubReg_Cellids.(subreg_names{k}).Limit; % on the limit of subregion
             
             
             if ismember(cid,csin);
                 if ~isnan(props(ic)); % if cell has already been assigned to one another subregion
                     fprintf('%s multiple loc: %d and %d...',cellids{ic},props(ic),k);
                     if rem(props(ic),1)~=0 % if prev. assigned subreg. was only on the Limit of sugreg, assign to this one
                         props(ic) = k;
                     else % if prev. assigned subreg. was IN the subreg., assign to the one with close centroid
                         ds = [props(ic) k];
                         [~,kk] = min(dist(ds));
                         props(ic) = ds(kk);
                     end;
                     fprintf('...at the end: %d\n',props(ic));
                 else
                     props(ic) = k;
                 end
                 
                 
             elseif ismember(cid,cslim);
                 
                 if ~isnan(props(ic)); 
                     fprintf('%s multiple loc: %d and %d...',cellids{ic},props(ic),k+0.2);
                     if rem(props(ic),1)~=0  % if prev. assigned subreg. was only on the Limit of sugreg, assign to the one with close centroid
                         ds = [props(ic)-0.2 k];
                         [~,kk] = min(dist(ds));
                         props(ic) = ds(kk)+0.2;
                     end % if prev. assigned subreg. was IN the subreg, choose that one
                     fprintf('...at the end: %d\n',props(ic));
                 else
                     props(ic) = k+0.2;
                 end;
             end
             
         end
     end
    
    
elseif contains(p.prop,'latency')
    
    if length(p.cellids)==2 
        props{1} = get_peaklats(p.cellids{1},['activation_' p.prop(1:strfind(p.prop,'latency')-2)],p.Events2Align,'sigma',p.sigma,'test_win',p.test_win);
        props{2} = get_peaklats(p.cellids{2},['inhibition_' p.prop(1:strfind(p.prop,'latency')-2)],p.Events2Align,'sigma',p.sigma,'test_win',p.test_win);
    else
        for rt = 1:length(p.resptype)
            props{rt} = get_peaklats(p.cellids,[p.resptype{rt} '_' p.prop(1:strfind(p.prop,'latency')-2)],p.Events2Align,'sigma',p.sigma,'test_win',p.test_win);
        end
    end
elseif contains(p.prop,'amplitude') || contains(p.prop,'rate')
    if length(p.cellids)==2 % activated vs inhibited cells
        props{1} = get_peakamps(p.cellids{1},'maxvalue',p.Events2Align,'sigma',p.sigma,'test_win',p.test_win);
        props{2} = get_peakamps(p.cellids{2},'minvalue',p.Events2Align,'sigma',p.sigma,'test_win',p.test_win);
    else
        for rt = 1:length(p.resptype)
            props{rt} = get_peakamps(p.cellids,p.resptype{rt},p.Events2Align,'sigma',p.sigma,'test_win',p.test_win);
        end
    end
    
elseif contains(lower(p.prop),'resp')
    
    test_window = [0 1];
    props = [];
    for evi = 1:size(p.Events2Align,2)
        event = p.Events2Align{evi};
        
        % Find responsive cells
        propname_resp = [event 'psth_stat_' num2str(test_window)];
        px = getvalue(propname_resp,p.cellids);
        props = cat(2,props,px);

    end
    
elseif contains(lower(p.prop),'pred')
    
    statgr_dir = fullfile(cell_dir,'grouped2','psth_stat1');
    load(fullfile(statgr_dir,'PredCells.mat'))
    
    props = [];
    for evi = 1:size(p.Events2Align,2)
        
        event = p.Events2Align{evi};
        px = nan(length(p.cellids),1);
        
        for ci = 1:length(p.cellids)
            act_cellid = p.cellids{ci};
            if ismember(act_cellid,PredCells.([event '_StopPartition']).none.F_smaller_than_S) % F<S
                px(ci) = 1;
            elseif ismember(act_cellid,PredCells.([event '_StopPartition']).none.F_bigger_than_S) % F>S
                px(ci) = -1;
            else
                px(ci) = 0;
            end
        end
        
    end
        props = cat(2,props,px);
        
elseif  contains(lower(p.prop),'patientgroup')
    
    if contains(lower(p.prop),'clinical');
        [patgroups, groups_nm] = clinical_groups;
    elseif contains(lower(p.prop),'RTchange');
        [patgroups, groups_nm] = clinical_groups({'RTdecrease','RTincrease'});
    end
    
    props = nan(length(cellids),1);
    
    for ic = 1:length(cellids)
        for k = 1:length(groups_nm)
            pats = patgroups{k};
            
            for j = 1:length(pats)
                patcells = findcell('rat',pats{j});
                if ismember(cellids{ic},patcells);
                    props(ic) = k;
                end
            end
            
        end
    end
    
    
    
    
    
end


%--------------------------------------------------------------------------
function peakamp = get_peakamps(cellids,value,EventTypes,varargin)


prs = inputParser;
addRequired(prs,'cellids',@iscell);
addRequired(prs,'value',@ischar);
addRequired(prs,'EventTypes',@iscell);
addParameter(prs,'sigma',0.01,@isnumeric);
addParameter(prs,'test_win',[0 1],@isvector);
parse(prs,cellids,value,EventTypes,varargin{:});
pr = prs.Results;

%%
global stat_dir


evnr = length(EventTypes);
peakamp = nan(length(cellids),evnr);
for ei = 1:evnr
    alignevent = EventTypes{ei};
    
    for ic = 1:length(cellids)
        
        act_cellid = cellids{ic};
        cellidF = regexprep(act_cellid,'\.','_');
        
        if pr.sigma==0.01
            respfile = dir(fullfile(stat_dir,alignevent,[cellidF '*.mat']));
            load(fullfile(respfile(1).folder,respfile(1).name)); close(gcf);
            
            stats1 = stats1{1,1};
        else
            
            partition = 'all';   % partition trials
            wn = [-3 3];   % full raster window in seconds
            dt = 0.001;   % raster resolution, in seconds
            bwin = [-2.5 -1];   % baseline window for MW-test
            twin = [0 1];
            
            
            [~, ~, ~, ~, ~, stats1] = ...
                ultimate_psth(act_cellid,'trial',alignevent,wn,...
                'dt',dt,'sigma',pr.sigma,'parts',partition,'isadaptive',0,...
                'maxtrialno',Inf,'baselinewin',bwin,'testwin',pr.test_win,'relative_threshold',0.1,'display',false); % calculate psth
        end
        
        
        peakamp(ic,ei) = ((stats1.(value)-stats1.baseline)*100)/stats1.baseline;
    end
    
end



%--------------------------------------------------------------------------
function peaklat = get_peaklats(cellids,resppeak,EventsType,varargin)

prs = inputParser;
addRequired(prs,'cellids',@iscell);
addRequired(prs,'resppeak',@ischar);
addRequired(prs,'EventsType',@iscell);
addParameter(prs,'sigma',0.01,@isnumeric);
addParameter(prs,'test_win',[0 1],@isvector);
parse(prs,cellids,resppeak,EventsType,varargin{:});
pr = prs.Results;



%%
global stat_dir cell_dir

evnr = length(EventsType);
peaklat = nan(length(cellids),evnr);
for ei = 1:evnr
    alignevent = EventsType{ei};
    
    for ic = 1:length(cellids)
        
        act_cellid = cellids{ic};
        cellidF = regexprep(act_cellid,'\.','_');
        
        if pr.sigma==0.01
            respfile = dir(fullfile(stat_dir,alignevent,[cellidF '*.mat']));
            
            load(fullfile(respfile(1).folder,respfile(1).name)); close(gcf);
            
            peaklat(ic,ei) =  stats1{1,1}.(resppeak);
        else
            
            
            partition = 'all';   % partition trials
            wn = [-3 3];   % full raster window in seconds
            dt = 0.001;   % raster resolution, in seconds
            bwin = [-2.5 -1];   % baseline window for MW-test
            twin = [0 1];
            
            
            [~, spsth, ~, ~, ~, stats1] = ...
                ultimate_psth(act_cellid,'trial',alignevent,wn,...
                'dt',dt,'sigma',pr.sigma,'parts',partition,'isadaptive',0,...
                'maxtrialno',Inf,'baselinewin',bwin,'testwin',twin,'relative_threshold',0.1,'display',false); % calculate psth
            pause(0.5); close(gcf);
            time = wn(1):1/1000:wn(2);
            wininx = dsearchn(time',pr.test_win');
            
            % Find peaks
            if contains(resppeak,'activation')
                [pk,peakinx] = max(spsth(wininx(1):wininx(2)));
            elseif contains(resppeak,'inhibition')
                [pk,peakinx] = min(spsth(wininx(1):wininx(2)));
            end
            
            if contains(resppeak,'peak')
                
                %                 resdir = fullfile(cell_dir,'Peak',['spsth' num2str(pr.sigma) '_twin' num2str(pr.test_win)],resppeak);
                %                 if ~isfolder(resdir); mkdir(resdir); end;
                
                
                
                %                 fig = figure;
                %                 plot(time,spsth); hold on;                %
                %                 scatter(time(wininx(1)+peakinx),pk,100,'g','filled');
                %
                %                 saveas(fig,fullfile(resdir,[cellidF '.jpg']))
                %                 saveas(fig,fullfile(resdir,[cellidF '.fig']))
                %                 close(fig)
                
                
                peaklat(ic,ei) =  time(wininx(1)+peakinx);
                
            elseif contains(resppeak,'slope')
                
                
                ma_wn = 50;
                coeff = ones(1, ma_wn)/ma_wn;
                spsthma = filter(coeff, 1, spsth);
                fDelay = (ma_wn-1)/2;
                
                
                slp = spsthma([wininx(1):wininx(1)+peakinx]+round(fDelay));
                slp2 = spsth([wininx(1):wininx(1)+peakinx]);
                
                smtime =  time([wininx(1):wininx(1)+peakinx]);
                
                dy = diff(slp);
                dy2 = diff(slp2);
                if contains(resppeak,'activation')
                    [mdif,mix] = max(dy);
                elseif contains(resppeak,'inhibition')
                    [mdif,mix] = min(dy);
                end
                peaklat(ic,ei) =  time(wininx(1)+mix);
                %
                %                 resdir = fullfile(cell_dir,'Slope',['spsth' num2str(pr.sigma) '_twin' num2str(pr.test_win)],resppeak);
                %                 if ~isfolder(resdir); mkdir(resdir); end;
                %
                %                 fig = figure;
                %                 subplot(121)
                %                 plot(time,spsth); hold on;
                %                 plot(time(ma_wn:end),spsthma(ma_wn:end));
                %                 %
                %                 subplot(122);
                %                 yyaxis left
                %                 plot(smtime,slp); hold on;
                %                 yyaxis right
                
                %                 plot(smtime(2:end),dy); hold on; plot(smtime(2:end),dy2);
                %
                %                 scatter(time(wininx(1)+mix),mdif,100,'g','filled');
                %                 %                 pause(0.11); close(fig);
                %                 saveas(fig,fullfile(resdir,[cellidF '.jpg']))
                %                 saveas(fig,fullfile(resdir,[cellidF '.fig']))
                %                 close(fig)
                
            end
            
        end
        
        
        
        
        
    end
    
end



%--------------------------------------------------------------------------
function [MRLs Ps FTM]= get_PCvalues(cellids,pc_event,fr_band,PCrectype,winnr,downsamp,chanmean)




global cell_dir

if contains(PCrectype,'LFP')
    rectype = 'LFP';
    if chanmean
        chm = 'chanmean_allSTN';
    else
        chm = 'by-channel';
    end
    
    pcdir = fullfile(cell_dir,'PC',[rectype '_phase_coupling'],chm);
elseif contains(PCrectype,'EEG')
    rectype = 'EEG';
    chm = 'F4';
    
    pcdir = fullfile(cell_dir,'PC',[rectype '_phase_coupling'],chm);
end

if ~iscell(pc_event)
    event = pc_event; evty = pc_event;
elseif (iscell(pc_event) && length(pc_event)==1)
    event = pc_event{1}; evty = pc_event{1};
elseif  (iscell(pc_event) && length(pc_event)==2)
    event = pc_event{1}; evty = pc_event{2};
end

[MRLs Ps FTM] = deal([]);

plot_win = [-1.5 1.5];


cellidstit = cellfun(@(x) [x(1:14) '_' x(16:end)],cellids,'UniformOutput',false);
cellnr = length(cellids);


if strcmp(pc_event,'continu')
    for icc = 1:cellnr
        load(fullfile(pcdir,event,evty,[cellids{icc} '_PC_results.mat']))
        
        PCs.(cellidstit{icc}) = PC_results.Hilb_PC; PC_results = [];
    end
else
    load(fullfile(pcdir,event,evty,[fr_band '_' num2str(plot_win)],['PC_results_ds' downsamp '_' num2str(winnr) 'win.mat']))
    for icc = 1:cellnr
        try
            PCs.(cellidstit{icc}) = PC_results.Hilb_PC.(cellidstit{icc});
        catch
            fprintf('No cell %s\n',cellidstit{icc})
            continue
        end
    end
end

if exist('PCs')==0
    return;
end;
MRLs = nan(winnr,cellnr);
remcells = fieldnames(PCs);
for ci = 1:cellnr
    if ~ismember(cellidstit{ci},remcells)
        continue
    else
        celi = find(ismember(remcells,cellidstit{ci}));
    end
    if iscell(PCs.(remcells{celi}).ResVect)
        noe = cellfun(@(x) ~isempty(x),PCs.(remcells{celi}).ResVect);
        MRLs(noe,ci) = abs([PCs.(remcells{celi}).ResVect{noe}]);
    else
        MRLs(ci) = abs(PCs.(remcells{celi}).ResVect);
    end
end


Ps = nan(winnr,cellnr);
remcells = fieldnames(PCs);
for ci = 1:cellnr
    if ~ismember(cellidstit{ci},remcells)
        continue
    else
        celi = find(ismember(remcells,cellidstit{ci}));
    end
    if iscell(PCs.(remcells{celi}).Ray_P)
        noe = cellfun(@(x) ~isempty(x),PCs.(remcells{celi}).Ray_P);
        Ps(noe,ci) = [PCs.(remcells{celi}).Ray_P{noe}];
    else
        Ps(ci) = PCs.(remcells{celi}).Ray_P;
    end
end


FTM = nan(winnr,cellnr);
remcells = fieldnames(PCs);
for ci = 1:cellnr
    if ~ismember(cellidstit{ci},remcells)
        continue
    else
        celi = find(ismember(remcells,cellidstit{ci}));
    end
    if iscell(PCs.(remcells{celi}).ResVect)
        noe = cellfun(@(x) ~isempty(x),PCs.(remcells{celi}).ResVect);
        FTM(noe,ci) = angle([PCs.(remcells{celi}).ResVect{noe}]);
    else
        FTM(ci) = angle(PCs.(remcells{celi}).ResVect);
    end
end



%--------------------------------------------------------------------------
function  [cellids cellnrs stn_dist] = get_stndist(dist,SubRegList,STN_subregions,separateSTN)

SubReg_names = fieldnames(SubRegList);
subreg_nr = length(SubReg_names);

if contains(dist,'centroid')
    srnm = dist(1:strfind(dist,'centroid')-2);
end

% Get cellids
if separateSTN
    cellids = cat(2,SubRegList.(srnm).IN,SubRegList.(srnm).Limit);
else
    cellids = {}; cellnrs = nan(1,subreg_nr);
    for j = 1:length(SubReg_names)
        sm = SubReg_names{j};
        cellids = cat(2,cellids,SubRegList.(sm).IN,SubRegList.(sm).Limit);
        cellnrs(j) = numel(cat(2,SubRegList.(sm).IN,SubRegList.(sm).Limit));
    end
end

cellnr = length(cellids);

% Get distances corresponding to cellids
[ratname,session,tetrode,unit] = cellfun(@cellid2tags,cellids,'UniformOutput',false);
sides = arrayfun(@(x) fastif(strcmp(session{x}(end),'l'),'left','right'),1:cellnr,'UniformOutput',false);

if contains(dist,'centroid')
    sinx = find(ismember(srnm,SubReg_names));
    stn_dist = cell2mat(arrayfun(@(x) STN_subregions.Patients.(ratname{x}).(sides{x}).(['Ch' num2str(tetrode{x})]).Dist(sinx),1:cellnr,'UniformOutput',false));
    
elseif contains(dist,'vert')
    stn_dist = cell2mat(arrayfun(@(x) STN_subregions.Patients.(ratname{x}).(sides{x}).(['Ch' num2str(tetrode{x})]).vert,1:cellnr,'UniformOutput',false));
end

