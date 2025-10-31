function [wcoh_pats, wcavg, pats2,fff] = wcoh_avg(event,varargin)
% WCOH_AVG Average across patients
% WCOH_AVG(event,...) calculates and draws average MSWC map across patients. 

% Johanna Petra Szabó, 05.2025
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu
%%
prs = inputParser;
addRequired(prs,'event',@ischar)
addParameter(prs,'subevent','',@ischar)
addParameter(prs,'patients','all',@(x) ischar(x)|iscell(x)) % 'all' | 'tremor-dominant'|'akinetic-rigid'|'mixed' | 'RTdecrease'| 'RTincrease'
addParameter(prs,'region',{},@iscell) % 'all' | {'Motor'} | {'Associative'} | {'Limbic'}
addParameter(prs,'close2centr','all',@ischar) % 'all' | 'Motor' | 'Associative' | 'Limbic'
addParameter(prs,'maxpower','all',@ischar) % 'all' | 'beta' | 'high_delta'
addParameter(prs,'chanmean',1,@isnumeric)
addParameter(prs,'side','left',@(x) ischar(x)|iscell(x)) % 'bothside' | {'left','right'}
addParameter(prs,'iscell',true,@(x) islogical(x)| isnumeric(x)) % [] | true | false
addParameter(prs,'freqs',[1 100],@(x) isvector(x)|isnumeric(x));
addParameter(prs,'isfig',true,@islogical);
addParameter(prs,'downsamp',false,@islogical);
addParameter(prs,'alpha',NaN,@isnumeric);
addParameter(prs,'mcorrect','cluster',@ischar);
addParameter(prs,'plot_win',[-2 2],@isvector);
addParameter(prs,'baseline_win',[],@(x) isvector(x)|isnumeric(x));
addParameter(prs,'avgband',false,@islogical);
parse(prs,event,varargin{:})
g = prs.Results;

global rootdir figdir_pd

eegchan = 'F4';
time = -2 :1/250: 2; time = time(1:end-1);
pltime = g.plot_win(1):1/250:g.plot_win(2);
timinx = dsearchn(time',pltime');
if ~isempty(g.baseline_win)
    bastime = g.baseline_win(1):1/250:g.baseline_win(2);
    basinx = dsearchn(time',bastime');
    basinx = basinx- timinx(1)+1;
end

if isempty(g.subevent)
    figdir = fullfile(figdir_pd,'introp_EEG_LFP','WCoh',event);
else
    figdir = fullfile(figdir_pd,'introp_EEG_LFP','WCoh',[event '_' g.subevent]);
end

if ~isfolder(figdir); mkdir(figdir); end;

wcoh_structname = fullfile(figdir, ['WavCH_WIN' num2str(g.plot_win) '_BAS']);
% wcoh_structname = fullfile(figdir, ['WavCH_WIN' num2str(g.plot_win) '_BAS' num2str(g.baseline_win)]);
%
if g.downsamp
    load([wcoh_structname '_downsamp.mat']);
else
    load([wcoh_structname '.mat']);
end

load(fullfile(figdir,'WavCH_f.mat'));
load(fullfile(figdir,'WavCH_coi.mat'));


% wcoh_avg_structname = fullfile(figdir, ['WavCH_avg_WIN' num2str(g.plot_win) '_BAS' num2str(g.baseline_win) '.mat']);
% % wcoh_avg_structname = fullfile(figdir, ['WavCH_avg_WIN' num2str(g.plot_win) '_BAS.mat']);
%
% if exist(wcoh_avg_structname)~=2
%     WavCH_avg  = struct;
% else
%     load(wcoh_avg_structname);
% end


% Freqs
if ~isempty(g.freqs)
    f_inx = ff>g.freqs(1)&ff<g.freqs(2);
    fff = ff(f_inx);
else
    f_inx = 1:length(ff);
    fff = ff;
end

if strcmp(g.patients,'all')
    pats = fieldnames(WavCH);
else
    
    goodpats = fieldnames(WavCH);
    grpats = clinical_groups({g.patients},'intraop','stimoff');
    grpats = cat(2,grpats{:});
    pats = intersect(goodpats,grpats); % some patients might have missing data
end

if contains(g.side,'both')
    sids = {'left','right'};
else
    if ~iscell(g.side); sids = {g.side}; else; sids = g.side; end;
end

k = 1; wcoh = {}; pats2 = {};
for pp = 1:length(pats)
    patnm = pats{pp};
    for ss = 1:length(sids)
        sidenm = sids{ss};
        
        if ~isfield(WavCH.(patnm),sidenm); 
            fprintf('No %s side data for %s\n',sidenm,patnm);
            continue;
            
        end
        %% Select channels to average
        inchans = {};
        subreg_tit = '';
        %         % Get channels in predefined subregion
        %         if ~isempty(g.region{1});
        %             inchans0 = cell(1,length(g.region));
        %             for rr = 1:length(g.region)
        %                 try
        %                     inchans0{rr} = get_chan_subregion(patnm,sidenm,g.region{rr});
        %                 catch
        %                     disp('')
        %                 end
        %             end
        %             inchans = cat(2,inchans0{:});
        %             inchans1 = inchans;
        %             subreg_tit = [g.region{:}];
        %         end
        
        
        if ~g.chanmean
            % Get channel closest to centroid of one STN subreg.
            if  ~strcmp(g.close2centr,'all')
                load(fullfile(rootdir,'Channels_close2centroids.mat'));
                inchans = centr_tab{[patnm '_' sidenm(1)],g.close2centr};
                subreg_tit = ['Close to ' g.close2centr ' centr.'];
                
            end
            
            if  ~strcmp(g.maxpower,'all')
                load(fullfile(rootdir,'maxFREQ_channels.mat'));
                inchans =  maxFREQ_channels{[patnm '_' sidenm(1)],g.maxpower};
                subreg_tit = ['Max ' g.maxpower ' power'];
                
            end
            
            % Get channels with cells
            if ~isempty(g.iscell)
                cids = findcell('rat',patnm); % cellids for patient
                cellids = cids(contains(cids,sidenm(1))); % cellids on the current side of patient
                
                cell_chan = nan(1,length(cellids));
                for cc = 1:length(cellids)
                    [~,~,cell_chan(cc),~] = cellid2tags(cellids{cc});
                end
                inchans = arrayfun(@(x) ['Ch' num2str(x)],unique(cell_chan),'UniformOutput',0);
                
                if ~g.iscell
                    fdnm = fieldnames(WavCH.(patnm).(sidenm));
                    allchans = cellfun(@(x) x(strfind(x,'_')+1:end),fdnm,'UniformOutput',0);
                    
                    inchans00 = inchans;
                    inchans = allchans( ~ismember(allchans,inchans00) );
                end
                inchans2 = inchans;
                cell_tit = ['cells-' char(string(g.iscell))];
            else
                cell_tit = '';
            end
            
            % If no condition, get all channels
            if ~isempty(g.region{1})&& ~isempty(g.iscell)
                inchans = intersect(inchans1,inchans2);
                
                %         elseif isempty(g.region{1})&& isempty(g.iscell) && strcmp(g.close2centr,'all') && strcmp(g.maxpower,'all');
                %             fdnm = fieldnames(WavCH.(patnm).(sidenm));
                %             inchans = cellfun(@(x) x(strfind(x,'_')+1:end),fdnm,'UniformOutput',0);
            end
        else
            %             inchans = {'chanmean'};
            fdnm = fieldnames(WavCH.(patnm).(sidenm));
            inchans = cellfun(@(x) x(strfind(x,'_')+1:end),fdnm,'UniformOutput',0);
            
        end
        
        
        
        wcoh_chans = cell(1,length(inchans));
        for ci = 1:length(inchans)
            lfpchan = inchans{ci};
            if isfield(WavCH.(patnm).(sidenm),[eegchan '_' lfpchan])
                wcch = double(WavCH.(patnm).(sidenm).([eegchan '_' lfpchan]));
                
                if isempty(g.baseline_win)
                    % Remove outliers
                    wcch(abs(wcch)>1) = NaN;
                    wcch(isnan(wcch)) = nanmean(wcch(:));
                end
                
                wcoh_chans{ci} = wcch;
            else
                fprintf('No %s in %s %s\n',lfpchan,patnm,sidenm);
            end
        end
        if ~isempty(cat(3,wcoh_chans{:}))
            wcoh{k} =  mean(cat(3,wcoh_chans{:}),3);
            k = k+1;
            pats2 = cat(2,pats2,patnm);
        end
    end
end
wcoh_pats = cat(3,wcoh{:});

wcoh_pats = wcoh_pats(f_inx,timinx,:);

% 
% if ~isempty(g.baseline_win)
%     bas_avgx = repmat(mean(wcoh_pats(:,basinx,:),2), [1 size(wcoh_pats,2) 1]);
%     bas_stdx = repmat(std(wcoh_pats(:,basinx,:),[], 2,'omitnan'), [1 size(wcoh_pats,2 ) 1]);
%     wcoh_pats_indivnorm = (wcoh_pats - bas_avgx)./bas_stdx;
% else
%     wcoh_pats_indivnorm = wcoh_pats;
% end


if ~isempty(g.baseline_win)
    bas_avg = repmat(mean(wcoh_pats(:,basinx,:),[2 3]), [1 size(wcoh_pats,2) 1]);
    bas_std = repmat(std(wcoh_pats(:,basinx,:),[], [2 3],'omitnan'), [1 size(wcoh_pats,2) 1]);
    wcavg = (mean(wcoh_pats,3) - bas_avg)./bas_std;
% wcavg = mean(wcoh_pats_indivnorm,3);
else
    wcavg = mean(wcoh_pats,3);
end


nnr = size(wcoh_pats,3); % nr of patients included

% File name to save fig
if ~g.chanmean
    
    tit = {[g.patients ' patients (n=' num2str(nnr) ')'],event,['EEG: F4, LFP: ' subreg_tit ', ' cell_tit],[num2str(g.freqs) ' Hz']};
    fnm= [g.side '_AVG_' event '_' subreg_tit '_' cell_tit '_F' num2str(g.freqs) '_PATS_' g.patients '_DS' char(string(g.downsamp)) '_WIN' num2str(g.plot_win) '_BAS' num2str(g.baseline_win)];
else
    tit = {[g.patients ' patients (n=' num2str(nnr) ')'],event,['EEG: F4, LFP: chanmean'],[num2str(g.freqs) ' Hz']};
    fnm= [g.side '_AVG_' event '_chanmean_F' num2str(g.freqs) '_PATS_' g.patients '_DS' char(string(g.downsamp)) '_WIN' num2str(g.plot_win) '_BAS' num2str(g.baseline_win)];
end

if ~g.avgband
    if ~isnan(g.alpha)
%         data = permute(wcoh_pats_indivnorm,[3 1 2]);
%         zmap_cl_p = nan(size(data,[2 3]));
        
%         for fk = 1:length(fff)
%             [zmap_cl_p(fk,:), ~, ~,diffmap] = bas_permstat_pd(data(:,fk,:),1:length(basinx),500,g.alpha,g.alpha,g.mcorrect);
%         end
%         maskersp = zmap_cl_p~=0;
        
            [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(wcoh_pats,fff,g.alpha,1000,false,'fdr',[],1:length(basinx));
        
    end
    
    
    
    
    
    % WavCH_avg.(event).(['F' num2str(g.freqs)]).([subreg_tit '_' cell_tit]) = wcoh_avg;
    % save(fullfile(figdir,'WavCH_avg.mat'),'WavCH_avg');
    %%
    if g.isfig
        fig = figure;
        wcoh_fig(wcavg,[],pltime,fff,coicoi(timinx),[-.4 .4]);
        if ~isnan(g.alpha)
            hold on;
            zmap = norminv(exactp_ersp);
            h =  bootstatFDR_clustercorr(zmap,maskersp,pow2(log2(fff)),'intraop','EEG_LFP','Fp1',pltime,log2(fff))
%             contour(pltime, log2(fff),maskersp,'Color','red');
        end
        arrayfun(@(x) set(x,'Color','r'), h(~isnan(h)),'UniformOutput',0);
        title(tit);
        saveas(gcf,fullfile(figdir,[fnm '.jpg']))
        saveas(gcf,fullfile(figdir,[fnm '.fig']))
        %     saveas(gcf,fullfile(figdir,[fnm '.emf']))
        close(gcf);
    end
    
    
elseif g.avgband
    if g.isfig
        fig = figure;
        avgp = squeeze( mean( wcavg , 1) );
        sep = squeeze( (std( wcavg , [],1) ) / sqrt(size(wcoh_pats,1)));
        errorshade(pltime,avgp, sep, 'LineColor','k');
        ylim([-.6 .6])
        title(['Avg. across ' num2str(g.freqs) ' Hz'])
        ylabel('Norm. MSWC');
        xlabel(['Time rel. to ' event ' (sec)'])
        setmyplot_balazs(gca)
        title(tit);
        saveas(gcf,fullfile(figdir,[fnm '_bandAVG.jpg']))
        saveas(gcf,fullfile(figdir,[fnm '_bandAVG.fig']))
        saveas(gcf,fullfile(figdir,[fnm '_bandAVG.pdf']))
        close(gcf);
    end
    
end


end
