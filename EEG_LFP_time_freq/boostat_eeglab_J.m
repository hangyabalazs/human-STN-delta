function [exactp_ersp,maskersp,alpha2] = boostat_eeglab_J(P,f,alpha,naccu,isim,mcorrect,varargin)
%BOOTSTAT_EEGLAB_J      Bootstrap statistics on time-frequency data
%   [exactp_ersp,maskersp,alpha2] = boostat_eeglab_J(P,f,alpha,naccu,isim,mcorrect,...)
%       Performs bootstrap statistics using BOOTSTAT function of EEGLAB
%       toolbox (Delorme A & Makeig S (2004)).
%
%   Input parameters:
%       P               3D matrix: frequency x time x epochs
%
%       F               vector of frequency components
%
%       ALPHA           integer, alpha threshold to generate a mask map (1 where stat.
%                       sign., 0 otherwise)
%
%       NACCU           integer, number of exemplars to accumulate
%       
%       ISIM            true | false, if true time-frequency map with statistics plotted
%
%       MCORRECT        method to correct for multiple comparison 
%       
%   Optional parameters:
%       FORMULA         character array, formula to compute the given measure
%                       (default: 'mean(arg1,3);')
%
%       BASELINE_INX    time vector indices for baseline in second dimension, 
%                       if empty, all time points are used
%                       (default: [])
%
% See also: TIME_FREQ_PATIENTS

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu



if isempty(varargin)
    formula = 'mean(arg1,3);';
    basline_inx = [];
elseif nargin==7
    if isempty(varargin{1})
        formula = 'mean(arg1,3);';
    else
        formula = varargin{1};
    end
    basline_inx = [];
elseif nargin==8
    if isempty(varargin{1})
        formula = 'mean(arg1,3);';
    else
        formula = varargin{1};
    end
    basline_inx = varargin{2};
end

[ Pboot, ~, Pboottrials] = bootstat(P, formula, 'boottype', 'shuffle', ...
    'label', 'ERSP', 'bootside', 'both', 'naccu', naccu, ...
    'alpha', alpha/2,'dimaccu',[2],'basevect', basline_inx);


if isempty(basline_inx)
    PA = mean(P, 3);
else
    B = repmat(nanmean(P(:,basline_inx,:),[2 3]),[1 size(P,2)]);
    Bsd = repmat(std(P(:,basline_inx,:),[],[2 3],'omitnan'),[1 size(P,2)]);
    PA = (mean(P,3) - B) ./ Bsd; % Baseline correction for each trial for statistics
end

% pboot = permute( reshape(Pboottrials,[size(P,2) ceil(naccu/size(P,2)) size(Pboottrials,2)]), [3 1 2]);
Pboottrials = Pboottrials';

[exactp_ersp,~] = compute_pvals(mean(P,3), Pboottrials);
%%
if strcmp(mcorrect,'fdr')
    alpha2 = fdr(exactp_ersp, alpha);


else
    alpha2 = alpha;
end

maskersp = exactp_ersp <= alpha2;

if isim
    my_imagsc(PA,f,[]);
    hold on;
    contour(maskersp,'Color','white')
end

end


% reshaping data
% -----------
function [pvals,surrog] = compute_pvals(oridat, surrog, tail)

if nargin < 3
    tail = 'both';
end

if myndims(oridat) > 1
    if size(oridat,2) ~= size(surrog, 2) || myndims(surrog) == 2
        if size(oridat,1) == size(surrog, 1)
            surrog = repmat( reshape(surrog, [size(surrog,1) 1 size(surrog,2)]), [1 size(oridat,2) 1]);
        elseif size(oridat,2) == size(surrog, 1)
            surrog = repmat( reshape(surrog, [1 size(surrog,1) size(surrog,2)]), [size(oridat,1) 1 1]);
        else
            error('Permutation statistics array size error');
        end
    end
end

surrog = sort(surrog, myndims(surrog)); % sort last dimension

if myndims(surrog) == 1
    surrog(end+1) = oridat;
elseif myndims(surrog) == 2
    surrog(:,end+1) = oridat;
elseif myndims(surrog) == 3
    surrog(:,:,end+1) = oridat;
else
    surrog(:,:,:,end+1) = oridat;
end

[tmp idx] = sort( surrog, myndims(surrog) );
[tmp mx]  = max( idx,[], myndims(surrog));

len = size(surrog,  myndims(surrog) );
pvals = 1-(mx-0.5)/len;
if strcmpi(tail, 'both')
    pvals = min(pvals, 1-pvals);
    pvals = 2*pvals;
end;
end

function val = myndims(a)
if ndims(a) > 2
    val = ndims(a);
else
    if size(a,1) == 1,
        val = 2;
    elseif size(a,2) == 1,
        val = 1;
    else
        val = 2;
    end
end;
end