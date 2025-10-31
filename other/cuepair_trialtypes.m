function [ord all_ord ord_skip rev rev_skip] = cuepair_trialtypes(sessdir, TEfnm)
%CUEPAIR_TRIALTYPES sorts trials based on diffulty of presented cue pair.
% CUEPAIR_TRIALTYPES sorts trials into following groups based on cue pair:
% - ordered, i.e. numbers are in ascendent order: 12, 23 (level of conflict: lowest)
% - ordered but one number is skipped: 13 (level of conflict: middle)
% - reversed order, i.e. numbers are in descendent order: 32, 21 (level of conflict: middle)
% - reversed order and one number is skipped: 31 (level of conflict: highest)
%
% Input params:
%   TE: TrialEvents struct
%
% Johanna Petra Szabó, 10.2025
% Clinic for Neurosurgery and Neurointervention
% Semmelweis University, Budapest, Hungary
% szabo.johanna.petra@semmelweis.hu


TE = load(fullfile(sessdir,TEfnm));

if strcmp(fieldnames(TE),{'TE'})
    TE = TE.TE;
end


trnr = length( TE.StimulusOn );

cue1 = TE.Cue1; cue2 = TE.Cue2;
pairs = [cue1; cue2];

ord = []; ord_skip = []; rev = []; rev_skip = []; all_ord = []; rev_cp = []; ord_cp = [];
for k = 1:trnr
    if cue1(k)<cue2(k)
        if cue1(k)==1&&cue2(k)==2 % with this 2-3 pair is excluded...
            ord = [ord k]; % reading order
            ord_cp = [ord_cp pairs(:,k)];
        end
        if cue2(k)-cue1(k)==1
            all_ord = [all_ord k];
        elseif cue2(k)-cue1(k)==2
            ord_skip = [ord_skip k]; % reading order but skipped
        end
    elseif cue1(k)>cue2(k)
        if cue1(k)-cue2(k)==1
            rev = [rev k]; % reverse
        elseif cue1(k)-cue2(k)==2
            rev_skip = [rev_skip k];  % reversed and skipped
            
            rev_cp = [rev_cp pairs(:,k)];
        end
    end
end
