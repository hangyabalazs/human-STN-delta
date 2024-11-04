function PC_cell_level(plot_win,freqs, fr_names,PCdir,EventTypes,SubEventTypes,varargin)
%PC_CELL_LEVEL  Spike-phase coupling 
% PC_CELL_LEVEL(plot_win,freqs, fr_names,PCdir,EventTypes,SubEventTypes,...)
%   -Calculates mean vector length (MRL) and preferred phase of coupling
%   between spiking activity and LFP/EEG FREQS frequency oscillation (SPC).
%   SPC is calculated in time windows (PLOT_WIN) around a specific behavioral
%   events (EVENTTYPES) or subevents (SUBEVENTTYPES).
%
%   -MRL is calculated for phase values derived from by hilbert transformation
%   on band-pass filtered data, see SPIK_PHAS_EXTRACTION).
%   Phase values associated with spiking are plotted on polar histograms
%   Rayleigh's Z-Statistic is calculated and included in plots.
%   P values are written below polar histograms with red, if significant
%   (<0.05).
%   Path to save results and plots is defined in PCDIR.
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
%                   (relevant only for intraop LFP data)
%                   true | false (default value: true)
%   'subevs'        if true, PC is calculated for trials separated based on
%                   subevents deinfed in SUBEVENTTYPES
%                   true | false (default value: false)
%   'subregion'     character array or cell array, uses channel data
%                   derived from listed STN subregions
%                   (relevant only for intraop LFP data)
%          'all' | 'Motor' | 'Limbic' | 'Associative'  (default value: 'all')
%   'isplot'        1 | 0, if 1 figures are generated and saved (def.
%                   value: 1)
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
%
% See also: SPIK_PHAS_EXTRACTION

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
addParameter(prs,'chanmean',true,@islogical);
addParameter(prs,'subevs',false,@islogical);
addParameter(prs,'subregion','all',@ischar);
addParameter(prs,'isplot',true,@islogical);
addParameter(prs,'downsamp','no',@ischar);
addParameter(prs,'dominantfreq',true,@islogical);
addParameter(prs,'PC_win',1,@isvector);
addParameter(prs,'cellids',{},@iscell);
addParameter(prs,'rectype','LFP',@ischar);
parse(prs,plot_win,freqs,fr_names,PCdir,EventTypes,SubEventTypes,varargin{:})
pr = prs.Results;




close all
[~,fold] = fileparts(PCdir);
frnr = length(fr_names);


if contains(fold,'LFP')&& pr.chanmean
    chtit = ['chanmean_' pr.subregion 'STN'];
elseif contains(fold,'LFP')&& ~pr.chanmean
    chtit = 'by-channel';
elseif strcmp(pr.rectype,'EEG')
    chtit = 'F4';
end


if pr.subevs; subs = [1 2]; else; subs = 1; end;


smwin = linspace(pr.plot_win(1),pr.plot_win(2),pr.PC_win+1);

if isempty(pr.cellids)
    allcellids = findcell;
else
    allcellids = pr.cellids;
end

%% Find minimum spikenr

if strcmp(pr.downsamp,'spike')
    MINspiknr = downsamp_spiknr(pr.EventTypes,pr.SubEventTypes,subs,pr.PC_win,smwin,pr.PCdir,chtit,allcellids);
end
%%


for fri = 1:frnr
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
            
            % Load extracted spike-coupled phase values
            resdir = fullfile(PCdir,chtit,event,evty);
            load(fullfile(resdir,'PC_phase_values_allcells.mat'));
            
            % Result directory for figures
            figdir = fullfile(resdir,[frnm '_' num2str(plot_win)],'Polar_hist_bycell');
            if ~isdir(figdir); mkdir(figdir); end;
            
            % Create/Complete (instead of overwriting) results mat
            if exist(fullfile(resdir,[frnm '_' num2str(plot_win)],['PC_results_ds' pr.downsamp '_' num2str(pr.PC_win) 'win.mat']))==2
                load(fullfile(resdir,[frnm '_' num2str(plot_win)],['PC_results_ds' pr.downsamp '_' num2str(pr.PC_win) 'win.mat']))
            else
                PC_results = struct;
            end
            
            
            
            mycellids = fieldnames(PC_phase_values_allcells);
            
            [Ray_P_value,Wat_P_uplim,Res_vect,MU,Best_fit] = deal(cell(length(mycellids),1));
            for ic = 1:length(mycellids)
                cellidF = mycellids{ic};
                
                fprintf('%s, %s...\n',cellidF,evty)
                
                cellid = cellidF;
                cellid(end-1) = '.';
                cellix = ismember(allcellids,cellid);
                
                if ~any(cellix);
                    continue;
                end
                
                phasvals =  PC_phase_values_allcells.(cellidF).(frnm).Hilb_phase_values;
                sp_ts0 = PC_phase_values_allcells.(cellidF).(frnm).spike_timestamps;
                
                if strcmp(pr.downsamp,'spike')
                    min_spiknr = MINspiknr(cellix);
                else
                    min_spiknr = 0;
                end
                
                if isnan(min_spiknr)
                    fprintf('Min spike is zero %s %s\n',cellidF,event);
                    continue
                end
                
                phas_smwins = cell(1,pr.PC_win);
                for w = 1:pr.PC_win
                    phas_smwin0 = phasvals(sp_ts0>=smwin(w)&sp_ts0<=smwin(w+1));
                    
                    if strcmp(pr.downsamp,'spike')
                        spiknr = length(phas_smwin0);
                        new_spinx = sort(randperm(spiknr,min_spiknr));
                        phas_smwins{w} = phas_smwin0(new_spinx);
                    else
                        phas_smwins{w} = phas_smwin0;
                    end
                end
                
                if any(cellfun(@length,phas_smwins)<5)
                    fprintf('%s spike nr<5\n',cellidF)
                    continue
                end
                
                fignm = [frnm '_' cellidF '_polar_ds_' pr.downsamp '_WIN' num2str(pr.PC_win)];
                
                
%                 [Ray_P_value{ic},Wat_P_uplim{ic},Res_vect{ic},MU{ic},Best_fit{ic}] =  ...
%                     PC_figstat3(phas_smwins,smwin,cellidF,evty,figdir,pr.isplot,[],fignm,'hist');
                
                [Ray_P_value{ic},Wat_P_uplim{ic},Res_vect{ic},MU{ic},Best_fit{ic}] =  ...
                    PC_figstat3(phas_smwins,smwin,cellidF,evty,figdir,pr.isplot,fignm,'hist');
                
                
                
                PC_results.Hilb_PC.(cellidF).Ray_P = Ray_P_value{ic};
                PC_results.Hilb_PC.(cellidF).Wat1_P = Wat_P_uplim{ic};
                PC_results.Hilb_PC.(cellidF).Mises_MU = MU{ic};
                PC_results.Hilb_PC.(cellidF).ResVect = Res_vect{ic};
                PC_results.Hilb_PC.(cellidF).Bestmodel = Best_fit{ic};
                
            end
            %             NonHomogen_Distr = cellfun(@(x) fastif(~isempty(x),x<=0.05,false),Ray_P_value);
            %             Normal_Distr = cellfun(@(x) fastif(~isempty(x),x<=0.05,false),Wat_P_uplim);
            % %             Best_fit = Best_fit';
            %             Tab = table(mycellids,NonHomogen_Distr,Normal_Distr);
            %
            %             PC_results.Table = Tab;
            
            save(fullfile(resdir,[frnm '_' num2str(plot_win)],['PC_results_ds' pr.downsamp '_' num2str(pr.PC_win) 'win.mat']),'PC_results');
        end
    end
end
end



%--------------------------------------------------------------------------
function MINspiknr = downsamp_spiknr(EventTypes,SubEventTypes,subs,PC_win,smwin,PCdir,chlab,allcellids)

subenr = length(subs);


spikNRs = nan(length(allcellids),length(EventTypes)*subenr,PC_win);


for ei = 1:length(EventTypes)
    event = EventTypes{ei};
    
    for sei = subs
        
        if length(subs)==1
            evty = event;
        else
            evty = SubEventTypes{ei,sei};
        end
        
        resdir = fullfile(PCdir,chlab,event,evty);
        enr = (ei-1)*subenr+sei;
        
        load(fullfile(resdir,'PC_phase_values_allcells.mat'));
        
        
        cellids = fieldnames(PC_phase_values_allcells);
        
        for ic = 1:length(cellids)
            cellidF = cellids{ic};
            cellid = cellidF;
            cellid(end-1) = '.';
            cellix = ismember(allcellids,cellid);
            
            frnms = fieldnames(PC_phase_values_allcells.(cellidF));
            frnm = frnms{1};
            sp_ts = PC_phase_values_allcells.(cellidF).(frnm).spike_timestamps;
            for w = 1:PC_win
                sps = sp_ts>=smwin(w)&sp_ts<smwin(w+1);
                spikNRs(cellix,enr,w) = sum(sps);
            end
        end
    end
end

MINspiknr = min(spikNRs,[],[2 3]);
MINspiknr(MINspiknr<15) = NaN;

end





