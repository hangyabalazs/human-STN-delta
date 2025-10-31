function all_resvect = PC_groups_f(plot_win,freqs,fr_names,PCdir,EventTypes,SubEventTypes,varargin)
%PC_GROUPS  Population phase and MRL distribution
% PC_GROUPS (plot_win,freqs,fr_names,PCdir,EventTypes,SubEventTypes,...)
%   plots population phase histogram (rose diagram) derived from preferred
%   coupling phases of a group of units (see name-value pair arguments).
%
% Required inputs:
%     PLOT_WIN          1x2 vector, time window relative to event timestamp in
%                       sec, to use for SPC calculation
%
%     FREQS             Nx2 matrix, boundaries of frequency bands of
%                       interest (N is number of frequency bands to analyse, first column is
%                       the lower limit, second column is the lower limit)
%
%     FR_NAMES          cell array, labels of frequency bands of interest
%                       (nr of cells has to correspond to nr of rows in FREQS)
%
%     PCDIR             path to save results and figures
%
%     EVENTTYPES        1xN cell array of event labels, ex: {'StimulusOn','StopSignal','KeyPress1','Feedback'};
%
%     SUBEVENTTYPES     Nx2 cell array of partition ("subevent") labels, each row
%                       corresponds to an event label, columns to partitions
%                       {'CueStim','StopStim';'FailedStopTrial','SuccesfulStopTrial';'CueResponse','StopResponse';'Correct','Error';};

% Optional input (name-value pairs with default values):
%   'chanmean'      if true, channel averaged LFP data is used, if false
%                   LFP data is plotted channel-by-channel
%                   (relevant only for intraop LFP data) (def. value:
%                   false)
%                   true | false (default value: true)
%   'subevs'        if true, PC is calculated for trials separated based on
%                   subevents deinfed in SUBEVENTTYPES
%                   true | false (default value: false)
%   'subregion'     character array or cell array, uses channel data
%                   derived from listed STN subregions
%                   (relevant only for intraop LFP data)
%          'all' | 'Motor' | 'Limbic' | 'Associative'  (default value: 'all')
%   'isplot'        true | false, if true figures are generated and saved (def.
%                   value: true)
%   'downsamp'      spike or trial nr. downsampled
%                    'no' | 'spike' (default value: 'no')
%   'dominantfreq'  if true, dominant frequency within predefined frequency
%                   band limits (FREQS) are used for PC calculation, see FIND_DOMINANT_FREQ_BANDS
%                   true | false (default value: true)
%   'PCwin'         nr of smaller time windows (time window of plot_win will
%                   be divided to PC_WIN nr of smaller windows to use for
%                   SPC calculation) (def. value: 1)
%   'LFPdownsamp'   numeric value (N), downsample LFP signal by keeping every
%                   N-th sample starting with the first (def. value: 10)
%   'cellids'       cell array of unit IDs to analyse, if empty all
%                   detected units are analyse (def. value: {}).
%   'rectype'       recording type to use for analyses ('EEG' | 'LFP')
%                   (def. value: 'LFP')
%   'group'         cell array of unit groups;
%        'all' | 'signPC' | 'signPC_LFP' ('signPC_LFP' relevant in case of EEG PD analysis, selects unit sign. coupled to LFP)
%       [EVENT ' resp'] | [EVENT ' pred'], (EVENT is a char. array of an event label)
%   'channelgroup'  cell array of channel groups; if not empty, selects units
%                   based on the channel they were detected on
%        '{}' | {'close2centr'} (def. value: {})
%   'side'          char. array, if not empty, only units detected at
%                   specified side are considered
%       'left' | 'right' | {} (def. value: {})
%   'components'    1:N double vector | [], nr. of components in von
%                   Mises distribution
%                   if not empty, von Mises distributions composed by
%                   the mixture of different number of components (COMPONENTS)
%                   are fitted on the population phase distribution (it uses an
%                   expectation maximalization algorithm and AIC, BIC,
%                   and parametric bootstrap p-value is calculated to find
%                   the best fitting model. (def. value: [1 2 3])
%
% See also: PC_FIGSTAT3, PC_CELL_LEVEL

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

prs = inputParser;
addRequired(prs,'plot_win',@ismatrix);
addRequired(prs,'freqs',@ismatrix);
addRequired(prs,'fr_names',@iscell);
addRequired(prs,'PCdir',@isdir);
addRequired(prs,'EventTypes',@iscell);
addRequired(prs,'SubEventTypes',@iscell);
addParameter(prs,'chanmean',false,@islogical);
addParameter(prs,'subevs',false,@islogical);
addParameter(prs,'subregion','all',@ischar);
addParameter(prs,'isplot',true,@islogical);
addParameter(prs,'downsamp','no',@ischar);
addParameter(prs,'dominantfreq',true,@islogical);
addParameter(prs,'PC_win',1,@isvector);
addParameter(prs,'rectype','LFP',@ischar);
addParameter(prs,'group','all',@ischar);
addParameter(prs,'channelgroup',{},@iscell);
addParameter(prs,'alpha',0.05,@isnumeric);
addParameter(prs,'components',[1 2 3],@(x) isvector(x)||isempty(x));
addParameter(prs,'side',{},@(x) ischar(x)|iscell(x));
parse(prs,plot_win,freqs,fr_names,PCdir,EventTypes,SubEventTypes,varargin{:})
pr = prs.Results;


close all;

[~,fold] = fileparts(PCdir);



if contains(fold,'LFP')&& pr.chanmean
    chtit = ['chanmean_' pr.subregion 'STN'];
elseif contains(fold,'LFP')&& ~pr.chanmean
    chtit = 'by-channel';
elseif strcmp(pr.rectype,'EEG')
    chtit = 'F4';
end


smwin = linspace(pr.plot_win(1),pr.plot_win(2),pr.PC_win+1);
if pr.subevs; subs = [1 2]; else; subs = 1; end;

global rootdir group_dir


%%

for fri = 1:length(fr_names)
    if pr.dominantfreq && ~contains(pr.fr_names{fri},'dom')
        frnm = ['dom_' pr.fr_names{fri}];
    else
        frnm = [pr.fr_names{fri}];
    end
    
    
    for ei = 1:length(pr.EventTypes)
        event = pr.EventTypes{ei};
        
        for sei = subs
            
            if ~pr.subevs
                evty = event;
            else
                evty = pr.SubEventTypes{ei,sei};
            end
            
            resdir = fullfile(PCdir,chtit,event,evty);
            load(fullfile(resdir,'PC_phase_values_allcells.mat'));
            
            load(fullfile(resdir,[frnm '_' num2str(plot_win)],['PC_results_ds' pr.downsamp '_' num2str(pr.PC_win) 'win.mat']));
            
            
            % Select cells in group
            if strcmp(pr.group,'all')
                cellids = fieldnames(PC_results.Hilb_PC);
                
                
            elseif contains(pr.group,'signPC')
                
                cellids = find_signPC_cellids(pr,PCdir,event,evty,plot_win,resdir,frnm,PC_results);
                
            elseif contains(pr.group,'resp')
                load(fullfile(group_dir,'RespCells.mat'));
                respevent =  pr.group(1:strfind(pr.group,'resp')-2);
                
                if ~strcmp(respevent,'all')
                    if contains(pr.group,'Activ')
                        cellids = RespCells.(respevent).none.Activ;
                    elseif contains(pr.group,'Inhib')
                        cellids = RespCells.(respevent).none.Inhib;
                    else
                        cellids = cat(2,RespCells.(respevent).none.Activ, RespCells.(respevent).none.Inhib);
                    end
                    
                    
                elseif strcmp(respevent,'all')
                    celli = cell(1,3);
                    for ri = 1:3
                        switch ri; case 1; respev = 'StimulusOn';
                            case 2; respev = 'StopSignal';
                            case 3; respev = 'KeyPress1';
                        end;
                        
                        
                        if contains(pr.group,'Activ')
                            celli{ri} = RespCells.(respev).none.Activ;
                        elseif contains(pr.group,'Inhib')
                            celli{ri} = RespCells.(respev).none.Inhib;
                        else
                            celli{ri} = cat(2,RespCells.(respev).none.Activ, RespCells.(respev).none.Inhib);
                        end
                    end
                    cellids = unique( cat(2, celli{:} ) );
                end
                
            elseif contains(pr.group,'pred')
                
                load(fullfile(group_dir,'PredCells.mat'));
                respevent =  pr.group(1:strfind(pr.group,'pred')-2);
                
                if contains(pr.group,'SuccesfulStop')
                    cellids = PredCells.([respevent '_StopPartition']).none.SuccesfulStop;
                elseif contains(pr.group,'FailedStop')
                    cellids = PredCells.([respevent '_StopPartition']).none.FailedStop;
                else
                    cellids = cat(2,PredCells.([respevent '_StopPartition']).none.SuccesfulStop, ...
                        PredCells.([respevent '_StopPartition']).none.FailedStop);
                end
                
            elseif contains(upper(pr.group),'UA')
                
                sign_cellids = find_signPC_cellids(pr,PCdir,event,evty,plot_win,resdir,frnm,PC_results);
                
                props = get_prop('SUA',sign_cellids);
                if strcmp(upper(pr.group),'SUA')
                    cellids = sign_cellids(props);
                elseif strcmp(upper(pr.group),'MUA')
                    cellids = sign_cellids(~props);
                end
            else
                [patgroups, groups_nm] = clinical_groups({pr.group},'intraop','stimoff','left');
                
                % Get cellids - patient groups
                cellids0 = {};
                pats = patgroups{1};
                for j = 1:length(pats)
                    patcells = findcell('rat',pats{j});
                    cellids0 = cat(2,cellids0,patcells);
                end
                
                
                cellids_signPC = find_signPC_cellids(pr,PCdir,event,evty,plot_win,resdir,frnm,PC_results);
                cellids = intersect(cellids0,cellids_signPC);
            end
            
            %% Select cells based on channel group
            if ~isempty(pr.channelgroup)
                
                
                cellids2 = cellids;
                cellids = cell(length(cellids2),1); nic = 1;
                for ic = 1:length(cellids2)
                    cid = cellids2{ic};
                    [patnm,sessnm,achan,~] = cellid2tags(cid);
                    
                    
                    if ~isempty(pr.side) && ~strcmp(sessnm(end),pr.side(1));
                        continue
                    end
                    lfp_act_chan = ['Ch' num2str(achan)];
                    
                    if strcmp('close2centr',pr.channelgroup{1}) && ~strcmp('all',pr.channelgroup{2})
                        load(fullfile(rootdir,'Channels_close2centroids.mat'));
                        mychan = centr_tab{[patnm '_' sessnm(end)],pr.channelgroup{2}};
                        
                        %                     elseif strcmp('maxpower',pr.channelgroup{1}) && ~strcmp('all',pr.channelgroup{2})
                        %                         load(fullfile(rootdir,'maxFREQ_channels.mat'));
                        %                         mychan =  maxFREQ_channels{[patnm '_' sessnm(end)],pr.channelgroup{2}};
                        
                    end
                    
                    if ismember(lfp_act_chan,mychan);
                        cellids{nic} = cid;
                        nic = nic+1;
                    end
                    
                end
                cellids(cellfun(@isempty,cellids)) = [];
                
                grnm = [pr.group '_' pr.channelgroup{1} '_' pr.channelgroup{2}];
            else
                grnm = [pr.group];
            end
            
            for ic = 1:length(cellids)
                cellids{ic}(end-1) = '_';
            end
            %% Get phase values
            
            [all_phasvals,  all_resvect] = get_resvect(cellids,PC_results,pr.PC_win);
            
            %% Polar plot + stat
            
            % Result directory
            
            structdir = fullfile(PCdir,chtit,event,evty,[frnm '_' num2str(plot_win)],...
                'Groups_modelselection');
            figdir = fullfile(structdir,grnm);
            if ~isdir(figdir); mkdir(figdir); end;
            
            if ~isempty(pr.side); grnm = [grnm '_' 0]; end;
            
            fignm = [ grnm  '_ds_' pr.downsamp '_win' num2str(pr.PC_win)];
            structnm = ['PC_groups_ds_' pr.downsamp '_win' num2str(pr.PC_win) '.mat'];
            
            if exist(fullfile(structdir,structnm))==2
                load(fullfile(structdir,structnm))
            else
                PC_groups = struct;
            end
            
            if contains(pr.group,'resp')
                
                [r_p,wat_pup,avg_resvec,avg_mu,bestmod,mufit,kappafit,pfit] = ...
                    PC_figstat3(all_resvect,smwin,[pr.group ' cells'],event,figdir,pr.isplot,fignm,'scatter',pr.components);
                
            else
                
                [r_p,wat_pup,avg_resvec,avg_mu,bestmod,mufit,kappafit,pfit] = ...
                    PC_figstat3(all_phasvals,smwin,[pr.group ' cells'],event,figdir,pr.isplot,fignm,'hist',pr.components);
                
            end
            grnm2 = grnm; grnm2(isspace(grnm)) = '_';
            PC_groups.(grnm2).cellids = cellids;
            PC_groups.(grnm2).phasevals = all_phasvals;
            PC_groups.(grnm2).resvect = all_resvect;
            PC_groups.(grnm2).avg_ResVect = avg_resvec;
            PC_groups.(grnm2).avg_Mises_MU = avg_mu;
            PC_groups.(grnm2).avg_Ray_P = r_p;
            PC_groups.(grnm2).avg_Wat1_P = wat_pup;
            PC_groups.(grnm2).avg_Bestmodel = bestmod;
            PC_groups.(grnm2).avg_mufit = mufit;
            PC_groups.(grnm2).avg_kappafit = kappafit;
            PC_groups.(grnm2).avg_pfit = pfit;
            
            
            save(fullfile(structdir,structnm),'PC_groups')
            
            
            
            
            
            
        end
    end
end
end




%--------------------------------------------------------------------------
function [all_phasvals  all_resvect] = get_resvect(cellids,PC_results,PC_win)

% Get central phase values
all_phasvals  = cell(1,PC_win);
all_resvect  = cell(1,PC_win);

for ic = 1:length(cellids)
    cellidF = cellids{ic};
    cellidF(strfind(cellidF,'.')) = '_';
    try
        res_vect = PC_results.Hilb_PC.(cellidF).ResVect;
    catch
        res_vect = nan; continue
    end
    for w = 1:PC_win
        if ic==1
            all_phasvals{w} =  nan(1,length(cellids));
            all_resvect{w} =  nan(1,length(cellids));
        end
        all_phasvals{w}(ic) = angle(res_vect(w));
        all_resvect{w}(ic) = res_vect(w);
    end
end
end


function cellids = find_signPC_cellids(pr,PCdir,event,evty,plot_win,resdir,frnm,PC_results)
if contains(pr.group,'LFP')
    if strcmp(pr.rectype,'LFP')
        fprintf('SignPC_LFP group only for intraop EEG\n');
        return
    end
    PCdir_lfp = PCdir;
    es = strfind(PCdir,'EEG');
    PCdir_lfp(es:es+2) = 'LFP';
    resdir_lfp = fullfile(PCdir_lfp,'by-channel',event,evty);
    load(fullfile(resdir_lfp,[frnm '_' num2str(plot_win)],['PC_results_ds' pr.downsamp '_' num2str(pr.PC_win) 'win.mat']));
end
cellids0 = fieldnames(PC_results.Hilb_PC);
if pr.PC_win==1
    ray_sign = structfun(@(x) x.Ray_P<=0.05,PC_results.Hilb_PC);
elseif pr.PC_win==2
    ray_sign1 = structfun(@(x) x.Ray_P(1)<=0.05,PC_results.Hilb_PC);
    ray_sign2 = structfun(@(x) x.Ray_P(2)<=0.05,PC_results.Hilb_PC);
    ray_sign = ray_sign1|ray_sign2;
end
cellids = cellids0(ray_sign);
for ic = 1:length(cellids)
    cellids{ic}(end-1) = '.';
end
if strcmp(pr.rectype,'EEG')
    load(fullfile(resdir,[frnm '_' num2str(plot_win)],['PC_results_ds' pr.downsamp '_' num2str(pr.PC_win) 'win.mat']));
    
end

end