function patient_groups_compare(EventTypes)
%PATIENT_GROUPS_COMPARE     Comparison of patient groups
% PATIENT_GROUPS_COMPARE(EventTypes) compares UPDRS scores, clinical
% parameters, postoperative stimulation location, bursting index, mean firing rate
% and spike-phase coupling measures in patients with reaction time increase
% vs. patients with reaction time decrease. 
%
% Input parameters:
%     EVENTTYPES        1xN cell array of event labels, ex: {'StimulusOn','StopSignal'};
%

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu


global group_dir
grtype = 'RTchange';  

if strcmp(grtype,'clingroups');
    [patgroups, groups_nm] = clinical_groups;
elseif strcmp(grtype,'RTchange');
    [patgroups, groups_nm] = clinical_groups({'RTdecrease','RTincrease'});
end

resdir = fullfile(group_dir,'Patgroups'); if ~isfolder(resdir);mkdir(resdir); end;
grnr =length(groups_nm);

% Get cellids - patient groups
cellids = cell(1,grnr);
for k = 1:grnr
    pats = patgroups{k};
    for j = 1:length(pats)
        patcells = findcell('rat',pats{j});
        cellids{k} = cat(2,cellids{k},patcells);
    end
end



% Compare UPDRS scores/ clinical parameters
updrs_patgroups_comp(groups_nm)





% Compare postop stimulation localization

if exist(fullfile(rootdir,'postop_stim_loc.mat'))~=2
    % copy variables manually from excel table
    fprintf('Create a cell array manually with patient codes (N number) with variable name: patients !')
    keyboard;
    fprintf('Create N x 3 matrix with X(med-lat), Y(ant-post), Z(dors-ventr) coordinates of stimulated electrodes with variable names: left, right !')
    keyboard;
    
    stimloc_xyz = table(left, right,'RowNames',patients);
    
    fprintf('Create N x 3 matrix with distances to Associative, Limbic and Motor centroids of stimulated electrodes with variable names: left, right !')
    keyboard;
    
    stimloc_centrdist = table(left, right,'RowNames',patients);
    
    
    postop_stim_loc.XYZ = stimloc_xyz;
    postop_stim_loc.centr_dist = stimloc_centrdist;
    postop_stim_loc.XYZ_labels = {'X_lat','Y_AP','Z_vert'};
    postop_stim_loc. centr_labels = {'Associative','Limbic','Motor'};
    
    save(fullfile(rootdir,'postop_stim_loc.mat'),'postop_stim_loc')
end



for k = 1:2
    switch k; case 1; value = 'centr'; case 2; value = 'XYZ'; end
    patgroups_stimloc(groups_nm,value)
end





% Compare bursting index/ ratio of bursting units
BIthr = 0.35; % threshold for bursting index
bursting_comp(cellids,groups_nm,grtype,BIthr,resdir)




% Compare firing rate of units
meanFR_comp(cellids,groups_nm,grtype,resdir)


% Compare responsive units
resp_comp(cellids,groups_nm,grtype,resdir,EventTypes,[.05 0.01])


% Compare spike-phase coupling
PC_comp(cellids,groups_nm,grtype,resdir,'all',EventTypes,'LFP',[0 15],[0 0.8],false)


end


%--------------------------------------------------------------------------
function meanFR_comp(cellids,groups_nm,grtype,resdir)

grnr =length(groups_nm);
meanFR = cell(1,grnr);
for k = 1:grnr
    meanFR{k} = getvalue('mean_FR',cellids{k});
end


mL = max(cellfun(@length,meanFR));
meanFR_c = cellfun(@(x) cat( 1,x,nan(mL-length(x),1) ),meanFR,'UniformOutput',0  );
meanFR_bp = cell2mat(meanFR_c);

fig = figure;
subplot(1,2,1)
boxplot(meanFR_bp,groups_nm)
ylabel('mean FR');

if grnr>2; combs = nchoosek(1:grnr,2); else; combs = [1 2]; end;

cnr = size(combs,1);
for c = 1:cnr
    [rsum_p(c),~] = ranksum(   meanFR_bp(:,combs(c,1)), meanFR_bp(:,combs(c,2))  )
end
subplot(1,2,2); axis off;

yL = ylim;
for c = 1:cnr
    if rsum_p(c)<=0.05;  col = 'r';  else;  col = 'k';  end;
    text(0, yL(end)-(yL(end)/cnr)*c,...
        [num2str(combs(c,1)) ' vs ' num2str(combs(c,2)) ': Ranksum-p = ' num2str(rsum_p(c))],'Color',col)
end


resdir2 = fullfile(resdir,'meanFR');
if ~isfolder(resdir2);mkdir(resdir2); end;
saveas(fig,fullfile(resdir2,['meanFR_' grtype '.jpg']));
saveas(fig,fullfile(resdir2,['meanFR_' grtype '.fig']));
close(fig)
end




%--------------------------------------------------------------------------
function bursting_comp(cellids,groups_nm,grtype,BIthr,resdir)

grnr =length(groups_nm);

burstinx = cell(1,grnr);
for k = 1:grnr
    burstinx{k} = getvalue('BursIndex',cellids{k});
end
mL = max(cellfun(@length,burstinx));
burstinx_c = cellfun(@(x) cat( 1,x,nan(mL-length(x),1) ),burstinx,'UniformOutput',0  );
burstinx_bp = cell2mat(burstinx_c);

fig = figure;
subplot(1,2,1)
boxplot(burstinx_bp,groups_nm);
set_my_boxplot(gca)
ylabel('Bursting index');



if grnr>2; combs = nchoosek(1:grnr,2); else; combs = [1 2]; end;

cnr = size(combs,1);
for c = 1:cnr
    [rsum_p(c),~] = ranksum(   burstinx_bp(:,combs(c,1)), burstinx_bp(:,combs(c,2))  );
end
subplot(1,2,2); axis off;

yL = ylim;
for c = 1:cnr
    if rsum_p(c)<=0.05;  col = 'r';  else;  col = 'k';  end;
    text(0, yL(end)-(yL(end)/cnr)*c,...
        [num2str(combs(c,1)) ' vs ' num2str(combs(c,2)) ': Ranksum-p = ' num2str(rsum_p(c))],'Color',col);
end



resdir2 = fullfile(resdir,'Bursting');
if ~isfolder(resdir2);mkdir(resdir2); end;
saveas(fig,fullfile(resdir2,['BI_' grtype '.jpg']));
saveas(fig,fullfile(resdir2,['BI_' grtype '.fig']));
saveas(fig,fullfile(resdir2,['BI_' grtype '.pdf']));
close(fig)

% Compare nr. of bursting units

burst_color = rgb('maroon');
noburst_color = rgb('gray');
colors = [burst_color; noburst_color];


Bnr = cellfun(@(x) sum(x>=BIthr), burstinx);
noBnr = cellfun(@(x) sum(x<BIthr), burstinx);



fig = figure;
stacked_bar_perc([Bnr' noBnr'],colors,'k')
legend({['BI >= ' num2str(BIthr)], ['BI < ' num2str(BIthr)]})
xticks(1:grnr); xticklabels(groups_nm);
ylabel('Nr. of units')


combs = nchoosek(1:grnr,2);
cnr = size(combs,1);
for k = 1:cnr
    X = table([Bnr(combs(k,1));noBnr(combs(k,1))],[Bnr(combs(k,2));noBnr(combs(k,2))],...
        'VariableNames',groups_nm(combs(k,:)),...
        'RowNames',{['BI >= ' num2str(BIthr)],['BI < ' num2str(BIthr)]});
    [~,pval,stats(k)] = fishertest(X);
    
    yL  = ylim; xL = xlim;
    if pval<=0.05;  col = 'r';  else;  col = 'k';  end;
    text(xL(end)*.6, yL(end)-k*12,...
        [groups_nm{combs(k,1)}(1:2) ' vs ' groups_nm{combs(k,2)}(1:2)  ': Fisher p = ' num2str(pval)],'Color',col);
end

saveas(fig,fullfile(resdir2,['BInr_stacked_' grtype '.jpg']));
saveas(fig,fullfile(resdir2,['BInr_stacked_' grtype '.fig']));
saveas(fig,fullfile(resdir2,['BInr_stacked_' grtype '.pdf']));
close(fig)


end




%--------------------------------------------------------------------------
function updrs_patgroups_comp(patgr_labels)

global rootdir figdir_pd

resdir = fullfile(figdir_pd,'Behav','Behav_groups'); if ~isfolder(resdir); mkdir(resdir); end;
load(fullfile(rootdir,'UPDRS_scores.mat'))

patgroups = clinical_groups(patgr_labels);

varnms = updrs_tab.Properties.VariableNames;

    grnr = length(patgr_labels);
    up = cell(grnr,1);
for k = 1:length(varnms)
    
    for j = 1:grnr
        up{j} = updrs_tab{patgroups{j},varnms{k}};
    end
    
    mL = max(cellfun(@length,patgroups));
    up2 = cellfun(@(x) [x;  nan(mL-length(x),1)] , up, 'UniformOutput',0);
    bp = cat(2,up2{:});
    
    fig = figure;
    boxplot(bp,patgr_labels)
    set_my_boxplot(gca)
    
    hold on;
    for b = 1:size(bp,2)
        scatter(ones(1,mL)*b+randi([-10 10],1,mL)*0.01,bp(:,b),[],[0 0 0],'filled'); hold on;
    end
    
    if grnr==2; combs = [1 2]; elseif grnr==3; combs = [1 2; 1 3; 2 3]; end;
    pp = [0.9, 0.8, 0.7];
    for kk = 1:size(combs,1)
        [pv,~,~] = ranksum( bp(:, combs(kk,1) ), bp(:, combs(kk,2) ) );
        if pv<0.05; col = 'r'; else; col = 'k'; end;
        yL = ylim;
        text(1,yL(end)*pp(kk),[upper(patgr_labels{combs(kk,1)}(1:5)) ' vs ' upper(patgr_labels{combs(kk,2)}(1:5)) ': p=' num2str(pv)],'Color',col);
        
    end
    
    
    
    ylabel(varnms{k})
    
    tt = cellfun(@(x)  x(1:5), patgr_labels,  'UniformOutput',0);
    ttt = [tt{:}];
    
    
    fnm  = fullfile(resdir, [varnms{k} '_' ttt]);
    saveas(fig,[fnm '.jpg'])
    saveas(fig,[fnm '.fig'])
    saveas(fig,[fnm '.pdf'])
    close(fig)
end
end




%--------------------------------------------------------------------------
function patgroups_stimloc(groups_nm,value)
global rootdir figdir_pd

resdir = fullfile(figdir_pd, 'Behav','Behav_groups');
if ~isdir(resdir); mkdir(resdir); end;

load(fullfile(rootdir,'postop_stim_loc.mat'))


patgroups = clinical_groups(groups_nm);
grnr = length(patgroups);

if contains(value,'XYZ')
    stimloc_tab = postop_stim_loc.XYZ;
    coord_labs = postop_stim_loc.XYZ_labels;
elseif contains(value,'centr')
    stimloc_tab = postop_stim_loc.centr_dist;
    coord_labs = postop_stim_loc.centr_labels;
end
pats = stimloc_tab.Properties.RowNames;

pinx = cellfun(@(x) find(ismember(pats,x)), patgroups, 'UniformOutput',0);

sides = {'left' ,'right'};
grlab = cellfun(@(x) x(1:5) , groups_nm, 'UniformOutput',0);
    grlab2  = [grlab{:}];
    
    fig = figure;
    for h = 1:2 % sides
        for j = 1:3 % coordinates
            subplot(2,3,(h-1)*3+j);
            
            coord = stimloc_tab.(sides{h})(:,j);
            
            
            coord_gr = arrayfun(@(x)  coord(pinx{x}) ,1:length(pinx),'UniformOutput',0);
            
            mL = max(cellfun(@length, coord_gr));
            bp = cellfun(@(x) [x; nan(mL-length(x), 1)] ,coord_gr, 'UniformOutput',0);
            bp = cat(2,bp{:});
            boxplot(bp,grlab);
            
            [pval, ~, ~] = ranksum(bp(:,1),bp(:,2));
            xL = xlim; yL = ylim;
            if pval<0.05; col = 'r'; else; col = 'k'; end;
            
            text(xL(1),yL(1)+diff(yL)*0.9,['p=' num2str(pval)], 'Color',col);
            
            ylabel(coord_labs{j});
            title(['postop. stim in ' sides{h} ' STN']);
            set_my_boxplot(gca);
            
        end
    end
    
    set(fig,'Position', get(0,'Screensize') );
    
    fnm = fullfile(resdir,['stimloc_' value '_' grlab2]);
    saveas(fig, [fnm '.jpg'] )
    saveas(fig, [fnm '.fig'] )
    saveas(fig, [fnm '.pdf'] )
    close(fig);


end




%--------------------------------------------------------------------------
function resp_comp(cellids,groups_nm,grtype,resdir,EventTypes_all,alphas)

grnr =length(groups_nm);

% Nr of all cells
cellnrs = cellfun(@length,cellids);

% Responsive cells
resptypes = {'activation','inhibition'};
test_window = [0 1];
wn = [-3 3];
saveplot_wn  = [-1.2 3];
for ei = 1:length(EventTypes_all)
    event = EventTypes_all{ei};
    
    
    
    switch event
        case 'StimulusOn'
            colors = [0 0 1; 0.3010, 0.7450, 0.9330];
            bwin = [-2.5 -1];
            psth_colors = winter(grnr);
        case 'StopSignal'
            colors = [0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250];
            bwin = [-3 -1.5];
            psth_colors = autumn(grnr);
        case 'KeyPress1'
            colors = [0.4940, 0.1840, 0.5560; 0.75, 0, 0.75];
            bwin = [-3 -1.5];
            psth_colors = cool(grnr);
        case 'Feedback'
            colors = [0 1 0; 0.4660, 0.6740, 0.1880 ];
            bwin = [1.5 3];
            psth_colors = summer(grnr);
    end
    
    piecolors = cat(1,colors,[0.5 0.5 0.5]);
    
    % Find responsive cells
    propname_resp = [event 'psth_stat_' num2str(test_window)];
    
    resp_cellids = cell(grnr,2);
    for k = 1:grnr
        resp_values = getvalue(propname_resp,cellids{k});
        
        cellidx1 = find(resp_values==1); % activated
        cellidx2 = find(resp_values==-1); % inhibited cells
        
        resp_cellids{k,1} = cellids{k}(cellidx1);
        resp_cellids{k,2} = cellids{k}(cellidx2);
    end
    resp_cellnrs = cellfun(@length,resp_cellids);
    
    % Res folder
    resdir_resp = fullfile(resdir,'Respcells',event);
    if ~isfolder(resdir_resp);mkdir(resdir_resp); end;
    
    
    % PSTH maps
    isfig = true;
    psth_n = cell(grnr,length(resptypes));
    for k = 1:grnr
        for respk = 1:length(resptypes)
            if isfig;  Ha = figure; end;
            
            [psth_n{k,respk}] = norm_psth_map1(resp_cellids{k,respk},resptypes{respk},event,...
                'baseline','indiv','basl_psth',{},'bwin',bwin, 'parts','all','cLim',[-6 6],...
                'bindex',[],'isfig',isfig,'sigma',0.06);
            if isfig
                
                xL = xlim;
                xlim(interp1(wn,xL,saveplot_wn));
                
                fnm = [grtype '_psthMAP_' groups_nm{k} '_' resptypes{respk}];
                saveas(Ha,fullfile(resdir_resp,[fnm '.jpg']))
                saveas(Ha,fullfile(resdir_resp,[fnm '.fig']))
                close(Ha);
            end
        end
    end
    
    % PSTH avg
    %     avgfig = figure;
    psth_colors_light = 1-(1-psth_colors)*0.5;
    
    if grnr>2; combs = nchoosek(1:grnr,2); else; combs = [1 2]; end;
    
    
    for c = 1:size(combs,1)
        %         subplot(2,grnr,grnr+c)
        
        avgfig = figure;
        p1 = cell(1,2); j = 1;
        for k = combs(c,:)
            
            subplot(2,1,1)
            [p1{j}, ~,txt1,~] = avg_psth_plot(psth_n{k,1},[], event, ...
                [psth_colors(k,:)]);
            
            xlim(saveplot_wn); ylim([-5 5]);
            
            subplot(2,1,2)
            [~, ~,~,txt2] = avg_psth_plot(psth_n{k,2},[], event, ...
                [ psth_colors_light(k,:)]);
            
            xlim(saveplot_wn); ylim([-5 5]);
            
            hold on; j = j+1;
            delete(txt1); delete(txt2)
        end
        if ~any(cellfun(@isempty,p1))
            legend([p1{1}(end),p1{2}(end)],groups_nm(combs(c,:)),'Location','southeast','Autoupdate','off')
        end
        
        % STAT
        yL = ylim; yL = flip(yL);
        time = wn(1):1/1000:wn(2);
        for rt = 1:2 % Activ/ Inhib
            statdat = {psth_n{combs(c,1),rt}', psth_n{combs(c,2),rt}'};
            
            for ai = 1:length(alphas)
                [~, pgroup, ~, ~, statsgroup, statsinter] = std_stat(statdat,'groupstats','on',...
                    'mode','fieldtrip','fieldtripmethod','montecarlo','mcorrect','cluster',...
                    'fieldtripalpha',alphas(ai),'fieldtripclusterparam',{'clusterstatistic','wcm'});
                subplot(2,1,rt)
                if ~isempty(pgroup)
                    sign_p = pgroup{1};
                    
                    
                    if ai==1
                        draw_signifpatch(time,sign_p,[0.8 0.8 0.8])
                        %             scatter(plottime(sign_p),repmat(yL(1)*0.9,1,sum(sign_p)),[],[0.3 0.3 0.3],'filled')
                    elseif ai==2
                        draw_signifpatch(time,sign_p,[0.5 0.5 0.5])
                    elseif ai==3
                        draw_signifpatch(time,sign_p,[0.3 0.3 0.3])
                    end
                end
                
            end
        end
        fnm = [grtype '_psthAVG_' groups_nm{combs(c,1)} '_' groups_nm{combs(c,2)}];
        saveas(avgfig,fullfile(resdir_resp,[fnm '.fig']))
        saveas(avgfig,fullfile(resdir_resp,[fnm '.jpg']))
        saveas(avgfig,fullfile(resdir_resp,[fnm '.pdf']))
        close(avgfig)
    end
    
    
    %%
    % Pies with cell nrs
    piefig = figure;
    for k = 1:grnr
        subplot(1,grnr,k);
        cellnr_pie([resp_cellnrs(k,1), resp_cellnrs(k,2) cellnrs(k)-sum(resp_cellnrs(k,:))],...
            piecolors,cat(1,resptypes(:),'non-resp'));
        title(groups_nm{k});
    end
    
    set(piefig,'Position',get(0,'Screensize'));
    fnm = [grtype '_pies'];
    saveas(piefig,fullfile(resdir_resp,[fnm '.fig']))
    saveas(piefig,fullfile(resdir_resp,[fnm '.jpg']))
    close(piefig)
    
    
end
end





%--------------------------------------------------------------------------
function PC_comp(cellids,groups_nm,grtype,resdir,cellgr,EventTypes_all,rectype,rlim_hist,rlim_mrl,chanmean,varargin)



if ~isempty(varargin)
    respevent = varargin{1};
    resptype = varargin{2};
else
    respevent ='';
    resptype = '';
end


grnr =length(groups_nm);
piecolors = [0.5 0 0.5; 0.5 0.5 0.5];
PCcolors = cool(grnr);
PClinecolors = bone(grnr+1);

%
if contains(cellgr,'resp')
    all_cellids = cellids;
    cellids = cell(1,grnr);
    
    % Find responsive cells
    test_window = [0 1];
    propname_resp = [respevent 'psth_stat_' num2str(test_window)];
    
    resp_cellids = cell(grnr,2);
    for k = 1:grnr
        resp_values = getvalue(propname_resp,all_cellids{k});
        
        if contains(resptype,'activ')
            cellidx = find(resp_values==1); % activated
        else
            cellidx = find(resp_values==-1); % inhibited cells
        end
        
        cellids{k} = all_cellids{k}(cellidx);
    end
    
    
    resdir_PC = fullfile(resdir,[rectype '_PC'],[cellgr '_' respevent '_' resptype]);
    if ~isfolder(resdir_PC);mkdir(resdir_PC); end;
    
    
else
    
    resdir_PC = fullfile(resdir,[rectype '_PC']);
    if ~isfolder(resdir_PC);mkdir(resdir_PC); end;
end


cellnrs = cellfun(@length,cellids);

if grnr>2; combs = nchoosek(1:grnr,2); else; combs = [1 2]; end;
cnr = size(combs,1);

for ei = 1:length(EventTypes_all)
    event = EventTypes_all{ei};
    
    
    
    
    mrls = cell(1,grnr);
    ftms = cell(1,grnr);
    for k = 1:grnr
        mrls{k} = get_prop('PC MRL sign',cellids{k},'Events2Align',{event},...
            'win',1,'rectype',rectype,'fr_band','dom_high_delta','chanmean',chanmean);
        
        ftms{k} = get_prop('PC FTM sign',cellids{k},'Events2Align',{event},...
            'win',1,'rectype',rectype,'fr_band','dom_high_delta','chanmean',chanmean);
    end
    sign_cids = cellfun(@(x) ~isnan(x),mrls, 'UniformOutput',0);
    signPC_cellnrs = cellfun(@sum ,sign_cids);
    
    
    
    
    %% Pies
    piefig = figure;
    for k = 1:grnr
        subplot(1,grnr,k);
        cellnr_pie([ signPC_cellnrs(k) cellnrs(k)-signPC_cellnrs(k)],...
            piecolors,{'sign PC','not PC'});
        title(groups_nm{k});
    end
    
    set(piefig,'Position',get(0,'Screensize'))
    fnm = [grtype '_' cellgr '_' event '_PIES'];
    saveas(piefig,fullfile(resdir_PC,[fnm '.jpg']))
    saveas(piefig,fullfile(resdir_PC,[fnm '.fig']))
    saveas(piefig,fullfile(resdir_PC,[fnm '.pdf']))
    close(piefig);
    
    %% Chi-square
    
    stackfig = figure;
    mytbl = cat(1,signPC_cellnrs,cellnrs-signPC_cellnrs)
    for c = 1:cnr
        x2 = [ones(cellnrs(combs(c,1)),1)*2; ones(cellnrs(combs(c,2)),1)*3];
        x1 = cat(1,sign_cids{combs(c,:)});
        
        [chi_tbl, chi2, chi_p, labels] = crosstab(x1,x2);
        
        subplot(1,cnr,c);
        bar(chi_tbl','stacked'); xticks(1:2); xticklabels(groups_nm(combs(c,:)));
        yL = ylim;
        if chi_p<.05; col = 'r'; else; col = 'k'; end;
        text(0.2,yL(2)*.9,['chi-square p=' num2str(chi_p)],'Color',col);
        legend({'no PC','sign PC'},'Location','northoutside');
    end
    fnm = [grtype '_' cellgr '_' event '_CHISQ'];
    saveas(stackfig,fullfile(resdir_PC,[fnm '.jpg']))
    saveas(stackfig,fullfile(resdir_PC,[fnm '.fig']))
    saveas(stackfig,fullfile(resdir_PC,[fnm '.pdf']))
    close(stackfig);
    
    
    %% Polar plots
    
    polarfig = figure;
    % Polar histogram
    
    [polim,ray_p,pop_mrl,pop_ftm] = deal(nan(1,grnr));
    for c = 1:cnr
        subplot(2,grnr+1,c)
        
        
        for k = combs(c,:)
            phs = ftms{k}(~isnan(ftms{k}))';
            if ~isempty(phs)
                [polim(k), sign, ray_p(k), pop_mrl(k), Z, pop_ftm(k)] = PC_polar(phs,20,'ray',true,false,PCcolors(k,:));
            else
                polim(k) = polarhistogram(phs,20)
            end
            hold on
        end
        rlim(rlim_hist)
    end
    
    
    
    % Polar scatter
    for c = 1:cnr
        subplot(2,grnr+1,grnr+c+1)
        
        
        for k = combs(c,:)
            phs = ftms{k}(~isnan(ftms{k}))';
            ms = mrls{k}(~isnan(mrls{k}))';
            posc = polarscatter(phs,ms,50,PCcolors(k,:),'filled');
            hold on;
            pl(k) = polarplot([angle(pop_ftm(k)) angle(pop_ftm(k))],[0 pop_mrl(k)],...
                'Color',PClinecolors(k,:),'LineWidth',5)
            rlim(rlim_mrl)
            
        end
        
    end
    
    
    
    sp = subplot(2,grnr+1,grnr+1);
    legend(sp,[polim pl], repmat(groups_nm,[1 2]),'Location','southoutside'); axis off;
    yL = ylim;
    for k = 1:grnr
        if ray_p(k)<0.05;  col = 'r';  else;  col = 'k';  end;
        text(0, yL(end)-(yL(end)/grnr)*k,[groups_nm{k} ': ray-p = ' num2str(ray_p(k))],'Color',col)
    end
    
    
    % Stat
    w2_p = nan(2,2); rsum_p = nan(2,1);
    for c = 1:cnr
        j = 1;
        [phsStat,msStat] = deal(cell(2,1));
        for k = combs(c,:)
            phsStat{j} = ftms{k}(~isnan(ftms{k}))';
            msStat{j} = mrls{k}(~isnan(mrls{k}))';
            j = j+1;
        end
        if ~any(cellfun(@(x) length(x)<18,phsStat))
            
            [u2, w2_p(c,:)] = b_watsontwo(phsStat{1}',phsStat{2}');
            [rsum_p(c), ~] = ranksum(msStat{1},msStat{2});
        else
            w2_p(c,:) = NaN; rsum_p(c) = NaN;
        end
    end
    sp = subplot(2,grnr+1,2*(grnr+1)); axis off;
    yL = ylim;
    for c = 1:cnr
        if w2_p(c,2)<=0.05;  col = 'r';  else;  col = 'k';  end;
        text(0, yL(end)-(yL(end)/cnr)*c,[num2str(c) ': W2-p = ' num2str(w2_p(c,2))],'Color',col)
    end
    
    set(polarfig,'Position',get(0,'Screensize'))
    fnm = [grtype '_' cellgr '_' event '_POLAR'];
    saveas(polarfig,fullfile(resdir_PC,[fnm '.jpg']))
    saveas(polarfig,fullfile(resdir_PC,[fnm '.fig']))
    saveas(polarfig,fullfile(resdir_PC,[fnm '.pdf']))
    close(polarfig);
    
    %% MRL boxplots
    
    
    
    mL = max(cellfun(@length,mrls));
    mrls_n = cellfun(@(x)  cat(1, x, nan(mL-length(x),1)),mrls,'UniformOutput',0);
    mrls_bp = cat(2,mrls_n{:});
    
    
    fig = figure;
    boxplot(mrls_bp,groups_nm);
    set_my_boxplot(gca)
    yL = ylim;
    
    for c = 1:cnr
        
        if rsum_p(c)<=0.05;  col = 'r';  else;  col = 'k';  end;
        text(0.5, yL(end)*(1-c/10),[num2str(combs(c,1)) ' vs ' num2str(combs(c,2))  ': Ranksum-p = ' num2str(rsum_p(c))],'Color',col)
    end
    
    fnm = [grtype '_' cellgr '_' event '_BOX'];
    saveas(fig,fullfile(resdir_PC,[fnm '.jpg']))
    saveas(fig,fullfile(resdir_PC,[fnm '.fig']))
    saveas(fig,fullfile(resdir_PC,[fnm '.pdf']))
    close(fig);
    
end

end
