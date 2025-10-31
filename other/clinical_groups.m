function  varargout = clinical_groups(varargin)
%CLINIAL_GROUPS sorts patients according to group labels
%   patgroups = CLINIAL_GROUPS                  sorts patients according to clinical groups; 
%   patgroups = CLINIAL_GROUPS(groups_nm)       sorts patients according to groups specified in GROUPS_NM; 
%   {patgroups,groups_nm}  = CLINIAL_GROUPS     sorts patients according to clinical groups and outputs corresponding group labels; 

% 
%       Input parameters (optional):
%           GROUPS_NM   1 x m cell array of group labels
% 
%       Output parameters (optional):
%           PATGROUPS   1 x m cell array, where m is the number
%                       of groups; each cell contains 1 x n cell array, where n is the
%                       number of patients included in that group
%           GROUPS_NM   1 x m cell array of group labels

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

if nargin<1
    groups_nm =  {'tremor-dominant','akinetic-rigid','mixed'};
elseif nargin<2
    groups_nm = varargin{1};
elseif nargin>1
    groups_nm = varargin{1};
    rectime = varargin{2};
    condition = varargin{3};
end

global rootdir figdir_pd

grnr = length(groups_nm);
patgroups = cell(1,grnr);

if ismember(groups_nm,{'RTdecrease','RTincrease'})
    load(fullfile(rootdir,'Behav_groups.mat'));
    
    for gri = 1:grnr
        patgroups{gri} = Behav_groups.(groups_nm{gri});
    end
    
% elseif contains(groups_nm,{'low_UPDRS','high_UPDRS'})
%     
%     load(fullfile(rootdir,'UPDRS_patgroups.mat'));
% 
%     for gri = 1:grnr
%         patgroups{gri} = UPDRS_patgroups.(groups_nm{gri});
%     end
    
elseif ismember(groups_nm,{'tremor-dominant','akinetic-rigid','mixed'})
    
    clin_table = readtable(fullfile(rootdir,'PD_SSRT_clinical_data.xlsx'),'Sheet','Sheet1');
    
    for gri = 1:grnr
        ginx = ismember(clin_table.PDSubtype,groups_nm{gri});
        patgroups{gri} = clin_table.CodeName(ginx);
    end
    
elseif ismember(groups_nm,{'RevSkip_slower', 'RevSkip_faster'});
    
        load(fullfile(rootdir,'Ordered_vs_ReversedSkipped_groups.mat'));

     for gri = 1:grnr
        patgroups{gri} = Ordered_vs_ReversedSkipped_groups.([rectime '_' condition]).(groups_nm{gri});
    end
    
end
varargout{1} = patgroups;
varargout{2} = groups_nm;