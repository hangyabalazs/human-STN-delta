function RT_allgroupin(allpat_cntr,conditions,allpats,whichRT)
%RT_ALLGROUPSIN1 draws boxplots of RT/SSDp0.5/SSRT values 
%   RT_ALLGROUPSIN1(allpat_cntr,conditions,allpats,whichRT) plots 
%   RT/SSDp0.5/SSRT values (ALLPAT_CNTR)of type WHICHRT, for each task
%   condition (CONDITIONS), including all patients included in ALLPATS.
%   Patients assigned to different clinical group/ RT change group are marked with
%   distinct colour/ marker.
% 
%   Input parameters:
%       ALLPAT_CNTR     n x 1 cell array, where n is the number
%                          of task conditions to analyse; each cell contains a p x 2 matrix,
%                          where p is the number of patients included; the first column
%                          corresponds to left side task, second column to right side task
%                          results
% 
%       CONDITIONS      1x2 cell array of task conditions to compare
%               -first column corresponds to recording time, second column to
%               DBS stimulation condition
%               -ex.:{'preop','stimoff';'intraop','stimoff';'postop','stimoff';'postop','stimon'};
% 
%       ALLPATS     cell array of patients
% 
%       WHICHRT     which type of RT to compare
%           'RT' | 'Go_RT' | 'FAlarm_RT' (see calc_RT.m)


% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

global figdir_pd

condnr = size(conditions,1);


sides = {'left','right'}; sdnr = length(sides);

rt_figdir = fullfile(figdir_pd, 'Behav',whichRT); if ~isdir(rt_figdir); mkdir(rt_figdir); end;

clingroups = clinical_groups({'tremor-dominant','akinetic-rigid','mixed'});
behavgroups = clinical_groups({'RTdecrease','RTincrease'});

clininx = cellfun(@(x) ismember(allpats,x) , clingroups,'UniformOutput', 0);
behavinx = cellfun(@(x) ismember(allpats,x) , behavgroups,'UniformOutput', 0);

behavcolors = {[0 0 1]; [1 0 0]};
clinshapes = {'diamond','square',"^"};

lab = arrayfun(@(x) [conditions{x,1} ' ' conditions{x,2}],1:condnr,'UniformOutput', 0);

fig = figure;

for sidi = 1:sdnr
    sid = sides{sidi};
    
    x_oneside = [];
    for coci = 1:condnr
        x_oneside = [x_oneside allpat_cntr{coci}(:,sidi)];
    end
    
    subplot(sdnr,1,sidi)
    boxplot(x_oneside,lab); hold on;
    set_my_boxplot(gca)
    
    mL = size(x_oneside,1);
    for k = 1:condnr
        rands = randi([-100 100],mL,1)*0.002;
        for b = 1:2
            bix = behavinx{b};
            s(b) = scatter(ones(sum(bix),1)*k + rands(bix) , x_oneside(bix,k) ,[],behavcolors{b},'filled' );
        end
        for c = 1:3
            cix = clininx{c};
            s(2+c) = scatter(ones(sum(cix),1)*k + rands(cix) , x_oneside(cix,k),100,'k', clinshapes{c} ,'LineWidth',1.5);
        end
    end
    legend(s,{'RTdecrease','RTincrease','tremor-dominant','akinetic-rigid','mixed'})
    if contains(lower(whichRT),'perf')
        ylabel([whichRT '(%)'])
    else
        ylabel([whichRT '(s)'])
    end
    setmyplot_balazs(gca)
end

set(fig,'Position',get(0,'Screensize'));
fnm = fullfile(rt_figdir,['All_groups_in1_' num2str(condnr) 'conds' ]);
saveas(fig,[fnm '.jpg'])
saveas(fig,[fnm '.fig'])
saveas(fig,[fnm '.pdf'])
close(fig)

end
