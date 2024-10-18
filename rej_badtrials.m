function [EEG_ep1, EEG_ep2] = rej_badtrials(EEG_ep1,EEG_ep2,curr_resdir)


global  ALLEEG CURRENTSET EEG

% Reject bad trials/data
%-------------------------

EEG = EEG_ep2;

if EEG.trials~=1
    if exist([curr_resdir filesep 'rejected_epochs.mat'])==2; ifrej = 1; else; ifrej = 0; end;
else
    if exist([curr_resdir filesep 'events_with_boundaries.mat'])==2; ifrej = 1; else; ifrej = 0; end;
end


if ifrej==0
    pop_eegplot(EEG,1,1,1,[],'srate',EEG.srate,'spacing',75,...
        'eloc_file',EEG.chanlocs, 'winlength',30,'dispchans',32,'events',EEG.event,...
        'plottitle', 'Reject bad trials/data');
    
    
    fig1 = gcf;
    input('Select trials/data to reject, if ready press REJECT button, than any write any character to command window.\n');
end


if EEG.trials~=1
    
    
    if ifrej==0
        allep = 1:size(EEG_ep1.epoch,2);
        rejected_eps = allep(~ismember(allep,[EEG.epoch.index]));
        save([curr_resdir filesep 'rejected_epochs.mat'],'rejected_eps');
    elseif ifrej==1
        load([curr_resdir filesep 'rejected_epochs.mat']);
        EEG = pop_rejepoch(EEG,rejected_eps);
    end
    
    EEG_ep1 = pop_rejepoch(EEG_ep1,rejected_eps);

    
    
else
    if ifrej==0
        EEG_ep1.event = EEG.event;
        events_with_boundaries = EEG.event;
        save(fullfile(curr_resdir, 'events_with_boundaries.mat'),'events_with_boundaries')
    elseif ifrej==1
        load([curr_resdir filesep 'events_with_boundaries.mat'])
        EEG_ep1.event = events_with_boundaries;
        EEG.event = events_with_boundaries;
    end
    
end
EEG_ep2 = EEG;