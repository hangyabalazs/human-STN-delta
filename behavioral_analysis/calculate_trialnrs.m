function calculate_trialnrs
% CALCULATED_TRIALNRS Trial nr.s for each patient

allpats = {'pd01','pd02','pd03','pd04','pd05','pd07','pd08','pd10','pd12','pd13','pd14','pd15','pd16'};

gotbl = nan(length(allpats),3);
stoptbl = nan(length(allpats),3);
for rT = 1:3
    switch rT
        case 1;
            rectype = 'LFP';
            rectime = 'intraop'; condi = 'stimoff'; 
        case 2;
            rectype  = 'EEG';  rectime = 'postop';
              condi = 'stimoff'; 
        case 3
            rectype  = 'EEG';  rectime = 'postop';
               condi = 'stimon'; 
    end
    
    
    sess2analyse = getdata2analyse(filesdir, 'rectype',rectype,...
        'rectime',rectime,'patients', 'allpatients', 'side','right', 'condition',condi);
    
    for si = 1:length(sess2analyse)
        patinx = ismember(allpats,sess2analyse(si).patient);
        try
        load(fullfile(sess2analyse(si).sessfolder,['TrialEvents_' condi '.mat']));
        catch
            fprintf('No TE data %d %d %d\n', sess2analyse(si).patient,rectime, condi);
        end
        
        gotrial = length(TE.StimulusOn);
        stoptrial = sum(TE.RealStopTrial,'omitnan');
        gotbl(patinx,rT) = gotrial;
        stoptbl(patinx,rT) = stoptrial;
    end
    
end