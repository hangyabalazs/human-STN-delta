function fish_p = unit_subregions(EventTypes)
%UNIT_SUBREGIONS  Localization of units within the STN
% UNIT_SUBREGIONS(EventTypes)
%   - creates a struct with localization data to be used for analysis, based
%   on an excel table, localized in rootdir \ 'STN_loc' \ : 
%       'stn_regiok_final_anonim.xlsx'   name of excel file containing STN localization data
%            'STN_tavolsagok'    	name of excel sheet containing microelectrode coordinates ('X_lat','Y_AP','Z_vert' axes) +
%                                   distance from STN subregion centroids ('Associative','Limbic','Motor')
%            'STN_regiok'           name of excel sheet containing within
%                                   which STN subregion the microelectrodes are located  ('Associative','Limbic','Motor')
%   - saves the struct with localization data as 'STN_loc.mat' in
%   directory: rootdir \ 'STN_loc' \
%   - 
%   - generates pie charts and stacked bar plots with counts of units localized in distinct STN subregions
%   (all recorded units pooled/ behaviorally relevant units pooled/ behaviorally relevant units separated)
%   - tests for localization bias of behaviorally responsive / predictive
%   unis in STN subregions (Fisher's exact test)
%   - correlates unit-presenting microelectrode coordinates and distances from subregion
%   centroids with bursting indeces of units
% 
% Input:
%     EVENTTYPES        1xN cell array of event labels, ex: {'StimulusOn','StopSignal','KeyPress1','Feedback'};


% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


global rootdir 

% Create & save localization mat files + list of cellids in each subregion
if exist(fullfile(rootdir,'STN_loc','STN_loc.mat'))~=2
    create_subregion_struct
end


%  Percentage of units localized in distinct STN subregions (pie chart)
all_cells_pie


%  Percentage of responsive/predictive units localized in distinct STN subregions (pie chart)
allresp_pie(EventTypes)

%  Percentage of responsive/predictive units localized in distinct STN subregions (pie chart + stacked bar plots)
respcells_loc(EventTypes,'resp')
respcells_loc(EventTypes,'pred')


% Test localization bias of behaviorally responsive/ predictive units in STN subregions
fish_p = nan(length(EventTypes), 6);
for ei = 1:length(EventTypes)
    event = EventTypes{ei};
    
    for ri = 1:6
        if ri>3 && ~strcmp(event,'StopSignal'); continue; end;
        switch ri;
            case 1;resptyp = {'Activ','Inhib'}; case 2; resptyp = {'Activ'}; case 3; resptyp = {'Inhib'};
            case 4; resptyp = {'F_smaller_than_S', 'F_bigger_than_S'}; case 5; resptyp = {'F_smaller_than_S'}; case 6; resptyp = {'F_bigger_than_S'};
                
        end;
        
        fish_p(ei,ri) =  respcells_loc_stat(event,resptyp,EventTypes);
    end
end



% Correlate localization coordinates/ distnace from subregion centroids with bursting index
bursting_correlations
end




%--------------------------------------------------------------------------
function create_subregion_struct

global rootdir
excelname = 'stn_regiok_final_anonim.xlsx'; % name of excel file containing STN localization data
distance_sheet = 'STN_tavolsagok'; % name of excel sheet containing microelectrode coordinates + distance from subnucleus centroids
subreg_sheet = 'STN_regiok'; % name of excel sheet containing the information: within which STN subregion the microelectrodes are located
excelpath = fullfile(rootdir,'STN_loc',excelname);


STN_loc = struct;

STN_loc.subreg_names = {'Associative','Limbic','Motor'};
STN_loc.coordinate_names = {'X_lat','Y_AP','Z_vert'};
opts = detectImportOptions(excelpath);
%% Read distance data
opts.Sheet = distance_sheet;
opts.VariableNames = {'PatientNr','Side','Channel_name','ChannelNr','lat','AP',...
    'vert','Associative-lat','Associative_AP','Associative_vert','Limbic_lat','Limbic_AP','Limbic_vert',...
    'Motor_lat','Motor_AP','Motor_vert'}; % content of excel
opts.DataRange = 'A1:P158';
T = readtable(fullfile(excelpath),opts);

patientnrs = sort(T.PatientNr(~isnan(T.PatientNr)));
patient_codenames = arrayfun(@(x) fastif(x<10,['pd0' num2str(x)],['pd' num2str(x)]),patientnrs,'UniformOutput',0);


for k = 3:length(T.PatientNr)
    if isnan(T.PatientNr(k))
        T.PatientNr(k) = T.PatientNr(k-1);
    end
end; clear k;

T.PatientCodeNr(1) = ""; T.PatientCodeNr(2) = "";
T.PatientCodeNr(3:end) = arrayfun(@(x) fastif(T.PatientNr(x)<10,['pd0' num2str(T.PatientNr(x))],['pd' num2str(T.PatientNr(x))]),3:length(T.PatientNr),'UniformOutput',0);


for k = 3:length(T.Side)
    if isempty(T.Side{k})
        T.Side{k} = T.Side{k-1};
    end
end; clear k;


%%
for i = 1:length(patientnrs)
    ptnm = patient_codenames{i};
    
    pinx = find(cellfun(@(x) ~isempty(x),strfind(T.PatientCodeNr,ptnm)));
    
    for s = 1:2
        switch s; case 1; side = 'left'; sinx = 1:6; case 2; side = 'right'; sinx = 7:12; end;
        
        p_sinx = pinx(sinx(1):sinx(end));
        
        STN_loc.Patients.(ptnm).(side).Subregions.Associative = [str2double(T.Associative_lat{p_sinx(1)}),str2double(T.Associative_AP{p_sinx(1)}),str2double(T.Associative_vert{p_sinx(1)})]; % 'X_lat','Y_AP','Z_vert'
        STN_loc.Patients.(ptnm).(side).Subregions.Limbic = [str2double(T.Limbic_lat{p_sinx(1)}),str2double(T.Limbic_AP{p_sinx(1)}),str2double(T.Limbic_vert{p_sinx(1)})]; % 'X_lat','Y_AP','Z_vert'
        STN_loc.Patients.(ptnm).(side).Subregions.Motor = [str2double(T.Motor_lat{p_sinx(1)}),str2double(T.Motor_AP{p_sinx(1)}),str2double(T.Motor_vert{p_sinx(1)})]; % 'X_lat','Y_AP','Z_vert'
        
        
        
        chans = T.ChannelNr(p_sinx(2:end));
        chnr = length(chans);
        
        for j = 1:chnr
            cinx = p_sinx(j+1);
            if isempty(chans(j)) || isnan(chans(j))
                continue
            else
                chnm = ['Ch' num2str(chans(j))];
                STN_loc.Patients.(ptnm).(side).Channels.(chnm).Coordinates = [str2double(T.lat{cinx}),str2double(T.AP{cinx}),str2double(T.vert{cinx})]; % 'X_lat','Y_AP','Z_vert'
                STN_loc.Patients.(ptnm).(side).Channels.(chnm).Distance_centroids = [str2double(T.Associative_lat{cinx}),str2double(T.Limbic_lat{cinx}),str2double(T.Motor_lat{cinx})]; % 'Associative','Limbic','Motor'
            end
        end
    end
end

%% Read subregion data
opts.Sheet = subreg_sheet;
opts.VariableNames = {'PatientNr','Side','Channel_name','ChannelNr','Associative','Limbic','Motor'}; % content of excel
opts.DataRange = 'A1:G131';
T = readtable(fullfile(excelpath),opts);

patientnrs = sort(T.PatientNr(~isnan(T.PatientNr)));
patient_codenames = arrayfun(@(x) fastif(x<10,['pd0' num2str(x)],['pd' num2str(x)]),patientnrs,'UniformOutput',0);


for k = 2:length(T.PatientNr)
    if isnan(T.PatientNr(k))
        T.PatientNr(k) = T.PatientNr(k-1);
    end
end; clear k;

T.PatientCodeNr(1) = ""; T.PatientCodeNr(2) = "";
T.PatientCodeNr(2:end) = arrayfun(@(x) fastif(T.PatientNr(x)<10,['pd0' num2str(T.PatientNr(x))],['pd' num2str(T.PatientNr(x))]),2:length(T.PatientNr),'UniformOutput',0);


for k = 2:length(T.Side)
    if isempty(T.Side{k})
        T.Side{k} = T.Side{k-1};
    end
end; clear k;


%%
for i = 1:length(patientnrs)
    ptnm = patient_codenames{i};
    
    pinx = find(cellfun(@(x) ~isempty(x),strfind(T.PatientCodeNr,ptnm)));
    
    for s = 1:2
        switch s; case 1; side = 'left'; sinx = 1:5; case 2; side = 'right'; sinx = 6:10; end;
        
        p_sinx = pinx(sinx(1):sinx(end));
        
        
        chans = T.ChannelNr(p_sinx);
        chnr = length(chans);
        
        for j = 1:chnr
            cinx = p_sinx(j);
            if isempty(chans(j)) || isnan(chans(j))
                continue
            else
                chnm = ['Ch' num2str(chans(j))];
                subreg_resps = {T.Associative{cinx} T.Limbic{cinx} T.Motor{cinx}}; % 'X_lat','Y_AP','Z_vert'
                for jk = 1:3
                    resp = subreg_resps{jk};
                    if contains(resp,'igen')
                        newr = 1;
                    elseif contains(resp,'nem')
                        newr = -1;
                    elseif contains(resp,'határ')
                        newr = 2;
                    elseif isempty(resp)
                        newr = 0;
                    end
                    STN_loc.Patients.(ptnm).(side).Channels.(chnm).Subregion(jk) = newr;
                end
            end
        end
    end
end


%% Save cell lists grouped by subregion


pdcells_tags = listtag('tetrode');
patnms = pdcells_tags(:,1);
sessnms = pdcells_tags(:,2);
chanms = pdcells_tags(:,3);


pats = fieldnames(STN_loc.Patients);
subregs = STN_loc.subreg_names;

% Add location data to CellBase Matrix


for sr = 1:length(subregs)
    propname{sr} = subregs{sr};
    if ~ismember(propname{sr},listtag('prop'))
        insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname{sr})
    end
end; clear sr;

for ti = 1:size(pdcells_tags,1)
    cellids = findcell('RAT',patnms{ti},'SESSION',sessnms{ti},'TETRODE',chanms{ti});
    
    switch sessnms{ti}(end)
        case 'l'
            sid = 'left';
        case 'r'
            sid = 'right';
    end
    reg = STN_loc.Patients.(patnms{ti}).(sid).Channels.(['Ch' chanms{ti}]).Subregion;
    
    for sr = 1:length(subregs)
        setvalue(cellids,propname{sr},reg(sr));
    end; clear sr;
    
end


% Make list of cells in every subregion


for si = 1:length(subregs)
    sreg = subregs{si};
    STN_loc.SubReg_Cellids.(sreg).IN = {};
    STN_loc.SubReg_Cellids.(sreg).Limit = {};
    
    allcells = findcell;
    
    nxIN = 1; nxLim = 1;
    for ci = 1:length(allcells)
        cellid = allcells{ci};
        val = getvalue(sreg,cellid);
        if val==1
            STN_loc.SubReg_Cellids.(sreg).IN{nxIN} = cellid;
            nxIN = nxIN+1;
        elseif val==2
            STN_loc.SubReg_Cellids.(sreg).Limit{nxLim} = cellid;
            nxLim = nxLim+1;
        end
    end
end

%% Save
save(fullfile(rootdir,'STN_loc','STN_loc.mat'),'STN_loc');
end






%--------------------------------------------------------------------------
function bursting_correlations

global rootdir cell_dir

resdir = fullfile(cell_dir,'STN_loc','Bursting_STN_loc');
if ~isdir(resdir); mkdir(resdir); end;

load(fullfile(rootdir,'STN_loc','STN_loc.mat'))
pdcells = findcell;
allSubReg_names = STN_loc.subreg_names;
coord_names = STN_loc.coordinate_names;
centr_pair = {'Associative','Motor'};

[patgroups, groups_nm] = clinical_groups({'RTdecrease','RTincrease'});

%%
cellnr = length(pdcells);
di =nan(2,cellnr);
[co,motor_co , assoc_co] =deal(nan(3,cellnr));
[grtag, sdtag] = deal(nan(1,cellnr));
for ic = 1:cellnr
    cellid = pdcells{ic};
    
    
    [ratname,session,tetrode,~] = cellid2tags(cellid);
    if strcmp(session(end),'l'); side = 'left'; else; side = 'right';end;
    
    % Distance from centroids
    dist = STN_loc.Patients.(ratname).(side).Channels.(['Ch' num2str(tetrode)]).Distance_centroids;
    
    six = ismember(allSubReg_names,centr_pair);
    centr_pair_nm = allSubReg_names(six);
    di(:,ic) = dist(six);
    
    % Channel coordinates
    coord = STN_loc.Patients.(ratname).(side).Channels.(['Ch' num2str(tetrode)]).Coordinates;
    co(:,ic) = coord;
    
    % Centroid coordinates
    cecoord = STN_loc.Patients.(ratname).(side).Subregions.Motor;
    motor_co(:,ic) = cecoord;
    
    cecoord = STN_loc.Patients.(ratname).(side).Subregions.Associative;
    assoc_co(:,ic) = cecoord;
    
    % Tag for patient group
    grtag(1,ic) = find(cellfun(@(x) ismember(ratname,x), patgroups));
    
    % Tag for side
    sdtag(1,ic) = strcmp(side,'left');
    
end
co2 = co;
%%
[props, ~] = get_prop('bursting',pdcells);


nopr = isnan(props);

% Centroid distance - Motor vs Associative + bursting
fig = figure;

rsc1 = rand([1,cellnr])/5;
rsc2 = rand([1,cellnr])/5;

c1 = di(1,:)+rsc1;
c2 = di(2,:)+rsc2;
scatter(c1(~nopr)',c2(~nopr)',80,'filled','CData',props(~nopr)); colormap(jet);
cb = colorbar; cl.Label.String = 'Bursting index';

hold on;
bs = scatter(c1(nopr),c2(nopr),80,'filled','k');

xlabel(['Distance from ' centr_pair_nm{1} 'centr.'])
ylabel(['Distance from ' centr_pair_nm{2} 'centr.'])
hold on;
ix1 = grtag==1;
ix2 = grtag==2;
ps1 = scatter(c1(ix1),c2(ix1),100,'ko');
ps2 = scatter(c1(ix2),c2(ix2),100,'kd');
legend([ps1 ps2 bs],[groups_nm, 'No BI calc.']);

set(gcf,'Position',get(0,'Screensize'));

fnm = fullfile(resdir,'A_M_centr_BI_scatter');
saveas(fig,[fnm '.jpg']);
saveas(fig,[fnm '.fig']);
saveas(fig,[fnm '.pdf']);
close(fig);


%    X = [ones(sum(~nopr),1) c1(~nopr)' c2(~nopr)'];
%             b = regress(props(~nopr),X);
%             [~,~,~,~,stats] = regress(props(~nopr),X);
%             pval = stats(3);
%             rho = stats(1);
%
%             f = fit( [c1(~nopr)' c2(~nopr)'], props(~nopr), 'poly23' );
% plot(f, [c1(~nopr)' c2(~nopr)'],  props(~nopr));

%%
%   idx = kmeans([x' props],2);
%
%    fig = figure;
%
% rsc1 = rand([1,cellnr])/5;
% rsc2 = rand([1,cellnr])/5;
%
% c1 = x(1,:)+rsc1;
% c2 = x(2,:)+rsc2;
%  scatter(c1(~nopr)',c2(~nopr)','filled','CData',props(~nopr)); colormap(jet); colorbar
%  hold on;
%   scatter(c1(nopr),c2(nopr),'filled','k');
%     scatter(c1(props<.35),c2(props<.35),'filled','k');
%
%  xlabel(centr_pair{1})
%  ylabel(centr_pair{2})
%  hold on;
% ix1 = idx==1;
% ix2 = idx==2;
% % ix3 = idx==3;
%   scatter(c1(ix1),c2(ix1),80,'ko');
%   scatter(c1(ix2),c2(ix2),80,'kd');
% %   scatter(c1(ix3),c2(ix3),80,'ks');

%% Correlations

% All units

for j = 1:3
    
    fig = figure;
    switch j;
        case 1;
            gpr = ~isnan(props);
            lab = 'all_units';
        case 2;
            gpr = props>=.35&(~isnan(props))
            lab = 'bursting_units_035_thr';
        case 3;
            halfdist = max(di(2,:))/2;
            gpr = (~isnan(props))&di(2,:)'<halfdist;
            lab = 'close2Mcentr_units_halfdist';
            
        case 4;
            halfdist = max(di(2,:))/2;
            gpr = props>=.35&(~isnan(props))&di(2,:)'<halfdist;
            lab = 'close2Mcentr_units_halfdist_035_thr';
    end
    
    subplot(1,2,1)
    polypredcicall_mod(di(1,gpr)',props(gpr),0.05,'robust');
    hold on;
    scatter(di(1,gpr)',props(gpr),50,'k','filled');
    xlabel('Dist. from centr.'); ylabel('bursting index');
    title([centr_pair_nm{1} ' centr.'])
    
    subplot(1,2,2)
    polypredcicall_mod(di(2,gpr)',props(gpr),0.05,'robust');
    hold on;
    scatter(di(2,gpr)',props(gpr),50,'k','filled');
    xlabel('Dist. from centr.'); ylabel('bursting index');
    title([centr_pair_nm{2} ' centr.'])
    
    
    tit = lab; tit(strfind(tit, '_')) = ' ';
    suptitle(tit)
    set(gcf,'Position',get(0,'Screensize'));
    
    fnm = fullfile(resdir,['corr_bursting_centrdist_' lab]);
    saveas(fig,[fnm '.jpg']);
    saveas(fig,[fnm '.fig']);
    saveas(fig,[fnm '.pdf']);
    close(fig);
    
end

%% Coordinates


co2(1,sdtag==0)= -co(1,sdtag==0); % revert lat. coordinate value of left sided units

for j = 1:4
    fig = figure;
    switch j;
        case 1;
            gpr = ~isnan(props);
            lab = 'all_units';
        case 2;
            gpr = props>=.35&(~isnan(props))
            lab = 'bursting_units_035_thr';
        case 3;
            halfdist = max(di(2,:))/2;
            gpr = (~isnan(props))&di(2,:)'<halfdist;
            lab = 'close2Mcentr_units_halfdist';
            
        case 4;
            halfdist = max(di(2,:))/2;
            gpr = props>=.35&(~isnan(props))&di(2,:)'<halfdist;
            lab = 'close2Mcentr_units_halfdist_035_thr';
    end
    
    for k = 1:3
        subplot(1,3,k)
        rsc = rand([1,cellnr])/5;
        
        cx = co2(k,:)+rsc;
        polypredcicall_mod(cx(gpr)',props(gpr),0.05,'robust');
        hold on;
        scatter(cx(gpr)',props(gpr),'filled','k');
        
        xlabel(coord_names{k})
        ylabel('bursting index')
        %     hold on;
        %     ix1 = grtag==1;
        %     ix2 = grtag==2;
        %     scatter(cx(ix1),props(ix1),80,'ko');
        %     scatter(cx(ix2),props(ix2),80,'kd');
        
        
    end
    tit = lab; tit(strfind(tit, '_')) = ' ';
    suptitle(tit)
    
    set(gcf,'Position',get(0,'Screensize'));
    
    fnm = fullfile(resdir,['corr_bursting_chancoords_' lab]);
    saveas(fig,[fnm '.jpg']);
    saveas(fig,[fnm '.fig']);
    saveas(fig,[fnm '.pdf']);
    close(fig);
    
    
end
end


%--------------------------------------------------------------------------
function all_cells_pie

global rootdir cell_dir

resdir = fullfile(cell_dir,'STN_loc','Cell_distrib'); if ~isfolder(resdir);mkdir(resdir); end;


load(fullfile(rootdir,'STN_loc','STN_loc.mat'));
subregs = STN_loc.subreg_names;
subnr = length(subregs);


cellids = findcell;
fig = figure;

locnrs = nan(1,subnr+1);
locnrs = find_nrs(cellids,subnr);



pie(locnrs);
labs = cat(2,subregs,'Outside of STN');
labs = arrayfun(@(x) cat(2,labs{x},[', n = ' num2str(locnrs(x))]),1:subnr+1,'UniformOutput',0);

legend(labs,'Location','southoutside');
title('All recorded cells')

fnm = fullfile(resdir,'all_cells_pie');
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);
end



%--------------------------------------------------------------------------
function allresp_pie(EventTypes)



global rootdir cell_dir group_dir

resdir = fullfile(cell_dir,'STN_loc','Cell_distrib'); if ~isfolder(resdir);mkdir(resdir); end;


load(fullfile(rootdir,'STN_loc','STN_loc.mat'));
subregs = STN_loc.subreg_names;
subnr = length(subregs);


allcellids = findcell;




load(fullfile(group_dir,'RespCells.mat'));
load(fullfile(group_dir,'PredCells.mat'));



pdcellids = cell(4,length(EventTypes));
for ei = 1:length(EventTypes)
    respevent = EventTypes{ei};
    
    
    pdcellids{1,ei} = RespCells.(respevent).none.Activ;
    pdcellids{2,ei} = RespCells.(respevent).none.Inhib;
    if strcmp(respevent,'StopSignal')
        pdcellids{3,ei} = PredCells.([respevent '_StopPartition']).none.F_smaller_than_S;
        pdcellids{4,ei} = PredCells.([respevent '_StopPartition']).none.F_bigger_than_S;
    end
    
end
cellids = unique(cat(2,pdcellids{:}));
rp_cnr = length(cellids);


% Resp vs non-resp pie
fig = figure;
pievec =  [rp_cnr length(allcellids)-rp_cnr];
pie(pievec);
labs = {'Resp/ Pred cells', 'No resp/pred cells'};
labs = arrayfun(@(x) cat(2,labs{x},[', n = ' num2str(pievec(x))]),1:2,'UniformOutput',0);

legend(labs)
fnm = fullfile(cell_dir,'resp_pred_pie');
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);



fig = figure;


[srs, ~] = get_prop('STN_loc',cellids);
locnrs = nan(1,subnr+1);
for r = 1:subnr
    locnrs(r)= sum(srs==r|srs==r+0.2);
end

locnrs(subnr+1) = length(cellids)-sum(locnrs(1:subnr));

pie(locnrs);
labs = cat(2,subregs,'Outside of STN');
labs = arrayfun(@(x) cat(2,labs{x},[', n = ' num2str(locnrs(x))]),1:subnr+1,'UniformOutput',0);

legend(labs,'Location','southoutside');
title('All resp/pred cells')

fnm = fullfile(resdir,'all_resppred_cells_pie');
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);



end




%--------------------------------------------------------------------------
function respcells_loc(EventTypes,type)

% bar plots + pie charts + saves table for stats

global rootdir group_dir cell_dir

if contains(type,'resp')
    load(fullfile(group_dir,'RespCells.mat'));
    resptypes = {'activation','inhibition'};
elseif contains(type,'pred')
    
    load(fullfile(group_dir,'PredCells.mat'));
    resptypes = {'F<S','F>S'};
end


resdir = fullfile(cell_dir,'STN_loc','Cell_distrib'); if ~isfolder(resdir);mkdir(resdir); end;
load(fullfile(rootdir,'STN_loc','STN_loc.mat'));
subregs = STN_loc.subreg_names;
subnr = length(subregs);

locnrs = zeros(2,subnr+1);
all_locnrs = []; all_labs = {};
for ei = 1:length(EventTypes)
    respevent = EventTypes{ei};
    
    if contains(type,'resp')
        pdcellids{1} = RespCells.(respevent).none.Activ;
        pdcellids{2} =RespCells.(respevent).none.Inhib;
    elseif contains(type,'pred')
        if ~strcmp(respevent, 'StopSignal'); continue; end;
        pdcellids{1} = PredCells.([respevent '_StopPartition']).none.F_smaller_than_S;
        pdcellids{2} = PredCells.([respevent '_StopPartition']).none.F_bigger_than_S;
    end
    
    
    % Pie
    fig = figure;
    for k = 1:2
        
        % Nr of resp units assigned to subregions
        locnrs(k,:) = find_nrs(pdcellids{k},subnr);
        
        
        subplot(1,2,k)
        % PIE
        pie(locnrs(k,:));
        
        % Legend label
        labs = cat(2,subregs,'Ambiguous loc.');
        labs = arrayfun(@(x) cat(2,labs{x},[', n = ' num2str(locnrs(k,x))]),1:subnr+1,'UniformOutput',0);
        legend(labs,'Location','southoutside');
        % Title
        title([respevent resptypes{k}])
    end
    saveas(fig,fullfile(resdir,[respevent '_' type '.jpg']))
    saveas(fig,fullfile(resdir,[respevent '_' type '.fig']))
    saveas(fig,fullfile(resdir,[respevent '_' type '.pdf']))
    close(fig);
    
    all_locnrs = cat(1,all_locnrs,locnrs);
    all_labs = cat(2,all_labs,{[respevent ' A'],[respevent ' I']});
end

fig = figure;

% STACKED BAR
bar(all_locnrs,'stacked')

% Labels and legend
xt = get(gca,'XTick');
set(gca,'XTick',xt,'XTickLabel',all_labs)
legend(cat(2,subregs,'Ambiguous loc.'))
set(fig,'Position',get(0,'Screensize'));
fnm = fullfile(resdir,['all_stacked_' type]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);

% PIE chart - all responsive/ predictive cells
fig = figure;
allsum = sum(all_locnrs,1);
pie(allsum)
labs = cat(2,subregs,'Ambiguous loc.');
labs2 = arrayfun(@(x) cat(2,labs{x},[', n = ' num2str(allsum(x))]),1:subnr+1,'UniformOutput',0);

legend(labs2)
fnm = fullfile(resdir,['all_pie_' type]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig);


% SAVE TABLE

% Total nr of units in each subreg
allcellls = findcell;
tot_locnrs = nan(1,subnr+1);
tot_locnrs = find_nrs(allcellls,subnr);

all_locnrs2 = cat(1,all_locnrs,tot_locnrs);

Unit_Labels = cat(2,all_labs,{'Total nr of units'})';


allx = cell(1,4);
for j = 1:4
    allx{j} = all_locnrs2(:,j); % assign unit nrs to subregions
end
locTAB =table(allx{:});

locTAB.Properties.RowNames = Unit_Labels';
locTAB.Properties.VariableNames = labs2;
save(fullfile(resdir,[type 'loc_TABLE.mat']),'locTAB')
end





%--------------------------------------------------------------------------
function fish_p = respcells_loc_stat(event,resptyp,EventTypes_all)


global rootdir cell_dir group_dir
resdir = fullfile(cell_dir,'STN_loc','Cell_distrib_stat',[event '_resp']); if ~isfolder(resdir);mkdir(resdir); end;
load(fullfile(rootdir,'STN_loc','STN_loc.mat'));
subregs = STN_loc.subreg_names;
subnr = length(subregs);


load(fullfile(group_dir,'RespCells.mat'));
load(fullfile(group_dir,'PredCells.mat'));
resptypess = {'Activ','Inhib'};
predtypess = {'F_smaller_than_S', 'F_bigger_than_S'};



% All resp cells
for ei = 1:length(EventTypes_all)
    respevent = EventTypes_all{ei};
    
    pdcellids{ei,1} = RespCells.(respevent).none.Activ;
    pdcellids{ei,2} = RespCells.(respevent).none.Inhib;
    
    if strcmp(respevent,'StopSignal')
        pdcellids{ei,3} = PredCells.([respevent '_StopPartition']).none.F_smaller_than_S;
        pdcellids{ei,4} = PredCells.([respevent '_StopPartition']).none.F_bigger_than_S;
    end
end

sevix = ismember(EventTypes_all,event); % index of event
ry = ismember(cat(2, resptypess,predtypess),resptyp); % index of event response type
ev_cellids = [pdcellids{sevix,ry}]; % cells responsive to selected event & resp type


for k = 1:2
    switch k
        case 1
            notevs = findcell;
            notev_label = 'All_other';
        case 2
            notevs = unique(cat(2,pdcellids{:}));
            notev_label = 'Other_resp';
    end
    
    [~,~,bi] = intersect(ev_cellids,notevs); % index of cells resp both to the selected event & other events
    noev_cellids = notevs;
    noev_cellids(bi) = []; % clear cells resp both to the selected event & other events from not-resp group
    
    [locs_ev0, ~] = get_prop('STN_loc',ev_cellids); % assigns one subregion to each unit where possible
    locs_ev = floor(locs_ev0);
    locs_ev(isnan(locs_ev)) = 0;
    locs_ev(locs_ev==2) = 0;
    
    isloc_ev = locs_ev~=0;
    
    
    notev_label2 = notev_label; notev_label2(strfind(notev_label,'_')) = ' ';
    
    
    [locs_nev0, ~] = get_prop('STN_loc',noev_cellids); % assigns one subregion to each unit where possible
    locs_nev = floor(locs_nev0);
    locs_nev(isnan(locs_nev)) = 0;
    locs_nev(locs_nev==2) = 0;
    
    %     isloc_nev = locs_nev~=0;
    
    %     %% Chi-square
    %     % Include all defined subregion
    %     x1 = [ ones(sum(isloc_ev),1)*4; ones(sum(isloc_nev),1)*5 ];
    %     x2 = [ locs_ev(isloc_ev); locs_nev(isloc_nev)  ];
    %
    %
    %
    %     [chi_tbl, chi2, chi_p, labels] = crosstab(x1,x2);
    
    %% Fisher exact
    
    locnrs = nan(2,subnr+1);
    locnrs(1,:) = find_nrs(ev_cellids,subnr); % event resp
    locnrs(2,:) = find_nrs(noev_cellids,subnr); % event no resp
    fish_tbl = locnrs(:,[1 3]);
    [~,fish_p,stats] = fishertest(fish_tbl);
    
    
    %% Assoc vs Motor
    fig = figure;
    bar(fish_tbl,'stacked');
    
    xticks(1:2); xticklabels({[event ' resp -' resptyp{:}], notev_label2})
    yL = ylim;
    %     if chi_p<0.05; col = 'r'; else; col = 'k'; end;
    %     text(0.2,yL(2)*.9,['chi-sq p=' num2str(chi_p) ],'Color',col);
    if fish_p<0.05; col = 'r'; else; col = 'k'; end;
    text(0.2,yL(2)*.8,['fishex p=' num2str(fish_p)], 'Color',col );
    legend(subregs([1 3]));
    
    fnm = fullfile(resdir,[event ' resp -' resptyp{:} '_vs_' notev_label '_Motor_vs_Assoc']);
    saveas(fig,[fnm '.jpg'])
    saveas(fig, [fnm '.fig'])
    saveas(fig, [fnm '.pdf'])
    close(fig);
    
    %     %% Chi-square
    %     % Include all defined subregion + ambiguous loc
    %     x1 = [ ones(length(ev_cellids),1)*4; ones(length(noev_cellids),1)*5 ];
    %     x2 = [ locs_ev; locs_nev  ];
    %     [chi_tbl, chi2, chi_p, labels] = crosstab(x1,x2);
    %
    %     fig = figure;
    %     bar(chi_tbl,'stacked');
    %     xticks(1:2); xticklabels({[event ' resp -' resptyp{:}], notev_label2})
    %     yL = ylim;
    %     if chi_p<0.05; col = 'r'; else; col = 'k'; end;
    %     text(0.2,yL(2)*.9,['chi-sq p=' num2str(chi_p) ],'Color',col);
    %     legend( cat(2,{'Ambiguous loc.'},subregs([1 3])) );
    %
    %     fnm = fullfile(resdir,[event ' resp -' resptyp{:} '_vs_' notev_label '_Motor_vs_Assoc_vs_Ambig']);
    %     saveas(fig,[fnm '.jpg'])
    %     saveas(fig,[fnm '.fig'])
    %     saveas(fig,[fnm '.pdf'])
    %     close(fig);
end
end


%--------------------------------------------------------------------------
function locnrs = find_nrs(cellids,subnr)

locnrs = nan(1,subnr);
[srs, ~] = get_prop('STN_loc',cellids); % assigns one subregion to each unit where possible
for r = 1:subnr
    locnrs(1,r)= sum(srs==r|srs==r+0.2); % rows: activ/inhib; columns: subregions
end
locnrs(1,subnr+1) = length(cellids)-sum(locnrs(1,1:subnr));



end