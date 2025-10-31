function h = bootstatFDR_clustercorr(zmap_p,mask,ff,rectime,rectype,act_chan,varargin)
% BOOTSTATFDR_CLUSTERCORR Cluster-based correction of stat. result obtained
% by testing significant change in time-frequency data relative to baseline.
%   bootstatFDR_clustercorr(zmap_p,mask,ff,rectime,rectype,act_chan,...)
%   loads presaved struct that contains cluster distribution under the null
%   (see: generate_clus_distr_TF.m). Defines cluster threshold (95th percentile 
%   of H0 clust. distr.), separately for different frequency bands. 
%   Finds clusters in ZMAP_P difference map obtained by testing relative to baseline. 
%   Removes clusters below threshold. Draws contours of significant
%   clusters.
% 
% See also: GENERATE_CLUS_DISTR_TF
%
% Johanna Petra Szabó, 10.2025
% Clinic for Neurosurgery and Neurointervention
% Semmelweis University, Budapest, Hungary
% szabo.johanna.petra@semmelweis.hu


if ~isempty(varargin)
    pltime = varargin{1};
    plfreq = varargin{2};
else
    pltime = 1:size(mask,2);
    plfreq = 1:length(ff);
end

global rootdir
load(fullfile(rootdir,'Cluster_distrib.mat'))

if strcmp(rectype,'EEG_LFP')
    rectit = 'EEG_LFP_wcoh';
else
    rectit =[rectime '_' rectype];
end
if strcmp(act_chan,'chanmean');
    freqs = cat(1,Cluster_distrib.(rectit).Ch1{:,1});
%    chs = fieldnames(Cluster_distrib.([rectime '_' rectype]));
%    CL_Thr = [];
%    for c = 1:length(chs)
%     CL_Thr = cat(1, CL_Thr,[Cluster_distrib.([rectime '_' rectype]).(chs{c}){:,3}]);
%    end
%    CL_Thr = mean(CL_Thr,1);
    CL_Thr = [Cluster_distrib.(rectit).(act_chan){:,3}];

else
    freqs = cat(1,Cluster_distrib.(rectit).(act_chan){:,1});
    CL_Thr = [Cluster_distrib.(rectit).(act_chan){:,3}];
end

h = nan(1,size(freqs,1));
for fr = 1:size(freqs,1)
%     f_ix =  find(ff>=(freqs(fr,1)-diff(freqs(fr,:))*.15)&ff<=(freqs(fr,2)+diff(freqs(fr,:))*.15));
       f_ix =  find(ff>=freqs(fr,1)&ff<=freqs(fr,2));
    if length(f_ix)<2; continue; end;
    
    % Linear indeces within freq band
    
    m2 = zmap_p;
    m2(f_ix,:) = 999;
    f_Linx = find(m2==999);
         
    
    cluster_thresh = CL_Thr(fr);
    islands = bwconncomp(mask);
    
    mask_cl = mask;
    for clusi = 1:islands.NumObjects
        %
        %         if numel(islands.PixelIdxList{clusi}) < cluster_thresh % if size of true cluster smaller then threshold, set to zero
        %             zmap_cl_p(islands.PixelIdxList{clusi}) = 0;
        %         end
        
        if sum(abs(zmap_p(islands.PixelIdxList{clusi}))) < cluster_thresh % if sum of true cluster smaller then threshold, set to zero
            mask_cl(islands.PixelIdxList{clusi}) = 0;
        elseif ~any(ismember(islands.PixelIdxList{clusi},f_Linx)) % get rid of clusters not located within the freq band
            mask_cl(islands.PixelIdxList{clusi}) = 0;
        end
        
        
    end
    
      hold on;
[~,h(fr)] = contour(pltime,plfreq,mask_cl,'Color','w');
end
    
    

%     % Look at clusters at the upper/ lower margins of the freq band
%     cnr = length(f_ix);
%     Lnr = numel(mask_clFR);
%     Frowix = 1:cnr:Lnr-cnr+1;
%     Lrowix = cnr:cnr:Lnr;
%     margCL = cellfun(@(x) any(ismember(x,[Frowix Lrowix])), islands.PixelIdxList ); % clusters touching the margin of the freq band
%     mask_cl2 = mask;
%     for kk = 1:islands.NumObjects
%         contour(1:size(mask_cl,2),f_ix,mask_clFR,'Color','w')
%     end
    
    
    
  
end