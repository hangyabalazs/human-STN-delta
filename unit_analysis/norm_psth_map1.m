function [psth_RR] = norm_psth_map1(cellids,resptype,alignevent, varargin);
%NORM_PSTH_MAP1 Normalized peri-event time histogram (PSTH) plot
% [psth_RR] = NORM_PSTH_MAP1(cellids,resptype,alignevent,...)
%   Calculates smoothed PSTHs for each unit in CELLIDS (for details see ULTIMATE_PSTH). 
%   Performs baseline normalization. Plots PSTH map and average.
%
% Required inputs:
%     CELLIDS       cell array of unit IDs
% 
%     RESPTYPE      type of units included in CELLIDS
%           'activation' | 'inhibition' | 'cluster'
% 
%     ALIGNEVENT    char. array of event label
% 
% Optional inputs (name-value pairs with default value):
%     'baseline'        type of baseline normalization
%           'indiv' | 'common' (default: 'indiv')
% 
%     'bwin'            1x2 vector of time limits of baseline window is sec
%                       (default value: [-3 -2])
% 
%     'basl_psth'       cell array, previously calculated PSTHs with same length, used
%                       for common baseline normalization ( (default value:[])
% 
%     'sigma'           smoothing kernel for the smoothed PSTH, in seconds (def. value: 0.04)
% 
%     'parts'           char. array, partition label (see PARTITION_TRIALS) (def. value: 'all')
% 
%     'cLim'            color axis limits for (def. value: [-50 50])
% 
%     'bindex'          property values corresponding to CELLIDS, for
%                       sorting, if empty sorting is based on max/ min of activation/
%                       inhibition response
% 
%     'grouptags'       cell array, if plotted CELLIDS are separated to
%                       groups, writes group labels on plots (def. value:{})
% 
%     'group_limit'     marks the group separation index (def. value:[])     
% 
%     'isfig'           true | false, if ture generates plots (def. value: true)
%
%     'parttags'        vector of tags associated with PARTS partition tag,
%                       defined in defineLabelsColors_pd.m (def. value: [1 2])
%
% Output parameter:
%   PSTH_RR      plotted PSTHs

% See also: ULTIMATE_PSTH, PARTITION_TRIALS, DEFINELABELSCOLORS_PD

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

prs = inputParser;
addRequired(prs,'cellids',@iscell)
addRequired(prs,'resptype',@ischar) % 'activation' | 'inhibition' | 'cluster'
addRequired(prs,'alignevent',@ischar)
addParameter(prs,'excl_beforeEv',[],@(x) ischar(x)|isnumeric(x)) % exclude spikes before this event
addParameter(prs,'baseline','indiv',@ischar) % 'indiv' | 'common'
addParameter(prs,'bwin',[-3 -2],@isvector)
addParameter(prs,'basl_psth',{},@iscell) % if baseline = 'common'
addParameter(prs,'sigma',0.04,@isnumeric)
addParameter(prs,'parts','all',@ischar)
addParameter(prs,'cLim',[-50 50],@isvector)
addParameter(prs,'bindex',[],@(x) isvector(x)||isnumeric(x))
addParameter(prs,'grouptags',{},@iscell) 
addParameter(prs,'group_limit',[],@isnumeric) % if groupstags ~isempty
addParameter(prs,'isfig',true,@islogical)
addParameter(prs,'parttags',[1 2],@(x) isvector(x)||isnumeric(x))
parse(prs,cellids,resptype,alignevent,varargin{:});
pr = prs.Results;





if isempty(pr.bindex);  bi = 1; else; bi = 0; end;


% resptype: activation | inhibition | cluster


% Time vector
wn = [-3 3]; % time window
dt = 0.001; % time resolution
time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector

btime = pr.bwin(1)*1000:dt*1000:pr.bwin(2)*1000;   % baseline time vector
[~, baseline_inx] = intersect(time,btime);   % baseline indices

if strcmp(pr.parts,'all')
    partnr = 1;
else
    partnr = 2;
    [mylabels, mycolors] = makeColorsLabels(@defineLabelsColors_pd,...
    {[pr.parts(2:end) '=' num2str(pr.parttags(1)) ],[pr.parts(2:end) '=' num2str(pr.parttags(2))]});

end

% Ultimate PSTH

R = cell(length(cellids),5);
for iC = 1:length(cellids)
    
    
    [R{iC,1} R{iC,2} R{iC,3} R{iC,4} R{iC,5}] = ultimate_psth(cellids{iC},'trial', alignevent,wn,...
        'dt',dt,'display',false,'sigma',pr.sigma,'parts',pr.parts,'isadaptive',0,...
        'maxtrialno',Inf,'baselinewin',pr.bwin,'testwin',[0 0.5],'relative_threshold',0.1,...
        'first_event',pr.excl_beforeEv);
    
    if partnr==2
        TE = loadcb(cellids{iC},'TrialEvents');   % load trial events
        [trials_all, tgs,~, ~, PartNum] = partition_trials(TE,pr.parts);
        
        if PartNum ==1&& contains(tgs{1},'2'); 
            trials_all = [{[]} trials_all(1)]; pn = 2; 
        elseif PartNum ==1&& contains(tgs{1},'1'); 
            trials_all = [trials_all(1) {[]}]; pn = 2; 
        end;
        trials = trials_all(pr.parttags);
        notr = cellfun(@length,trials)<2;
        
        if any(notr)
            R{iC,2}(notr,:) = nan(1,size(R{iC,2},2));
        end
    end
end

spsthR = R(:,2);
Rtags = R(:,4);



for spi = 1:partnr
    
    if partnr==2 && pr.isfig; subplot(2,1,spi);
    end;
    
    % Normalization
    if ~isempty(R)
        
        if partnr==1
            r = cell2mat(spsthR);
        else
            %             r = cell2mat(cellfun(@(x) x(pr.parttags(spi),:), spsthR,'UniformOutput',0));
            r = nan(size(spsthR,1),length(time));
            for ri = 1:size(spsthR,1) % loop over cells
                partinx = find(strcmp(Rtags{ri},[pr.parts(2:end) '=' num2str(pr.parttags(spi))]));
                if ~isempty(partinx)
                    r(ri,:) = spsthR{ri,:}(partinx,:);
                end
            end
        end
        
        if strcmp(pr.baseline,'indiv')
            
            bl = r(:,baseline_inx);   % baseline matrix
            repbas_mean = repmat( nanmean(bl,2),1,size(r,2) );
            repbas_std = repmat( nanstd(bl,[],2),1,size(r,2) );
            
        elseif strcmp(pr.baseline,'common')
            
            if isempty(pr.basl_psth)
                bl = r(:,baseline_inx);   % baseline matrix
            elseif ~isempty(pr.basl_psth)
                b = cat(1,pr.basl_psth{:});
                bl = b(:,baseline_inx);   % baseline matrix
            end
            
            
            repbas_mean =  repmat( nanmean(bl,[1 2]),size(r) );
            repbas_std = repmat( nanstd(bl,[],[1 2]),size(r) );
            
            
        end
        
        % PSTH with baseline normalization
        psth_R = ( r - repbas_mean ) ./ repbas_std;
            
            
            
        if partnr==2;
            psth_RR{spi} = psth_R;
        else
            psth_RR = psth_R;
        end
        
    else
        psth_R = [];
        psth_RR = [];
    end
    
    if pr.isfig
        if isempty(psth_R)
            return
        end
        cellidt = cellfun(@(x) regexprep(x,'\_','-'),cellids,'UniformOutput',0);
        
        if spi==1
            if isempty(pr.bindex)
                switch resptype
                    case 'activation'
                        [~, sinx] = sort(max(psth_R,[],2), 'descend');   % sort by maximum
                        
                    case 'inhibition'
                        [~, sinx] = sort(min(psth_R,[],2), 'ascend');   % sort by maximum
                        
                        
                    case 'cluster'
                        tw = [0 1000]; bw = [-2500 -1000];
                        twinx = dsearchn(time',tw'); bwinx = dsearchn(time',bw');
                        postw = median(psth_R(:,twinx(1):twinx(2)),2);
                        prew = median(psth_R(:,bwinx(1):bwinx(2)),2);
                        acts = find(postw>=prew);
                        inhs = find(postw<prew);
                        [~,sinx1] = sort(max(psth_R(acts,:),[],2), 'descend');
                        [~,sinx2] = sort(min(psth_R(inhs,:),[],2), 'ascend');
                        sinx = [acts(sinx1); inhs(sinx2)];
                    otherwise
                        sinx = 1:size(psth_R,1);
                end
            else
                [~, sinx] = sort(pr.bindex, 'descend');   % sort by maximum
                
            end
        end
        
        imagesc(psth_R(sinx,:))
        caxis(pr.cLim)
        
        if  ~isempty(pr.grouptags)
            xL = xlim;
            yL = ylim;
            line(xL,[pr.group_limit+0.5 pr.group_limit+0.5],'Color','r');
            if pr.group_limit~=0
                text(xL(1),yL(1)+0.5,pr.grouptags{1},'Color','w')
            end
            if pr.group_limit~= size(psth_R,1)
                text(xL(1),pr.group_limit+1,pr.grouptags{2},'Color','w')
            end
        end
        
        if partnr ==1
            title([alignevent ' ' resptype])
        else
            title(mylabels{spi})
        end
        
        colorbar
        yticks(1:length(sinx))
         yticklabels(cellidt(sinx))
        
        
        xticks(dsearchn(time',[-3000 -2000 -1000 0 1000 2000 3000]'))
        xticklabels(cellfun(@num2str,{-3000, -2000, -1000, 0, 1000, 2000, 3000},'UniformOutput',0))
        xlabel(['Time from ' alignevent '(ms)'])
        
    end
    
end

if pr.isfig
    set(gcf,'Position',get(0,'Screensize'));
end

