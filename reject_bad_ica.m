function [EEG_ep1] = reject_bad_ica(EEG_ep1,EEG_ep2, curr_resdir)

% Reject bad trials manually 
% Perform independent component analysis on high-pass filtered EEG data
% (EEG_ep2), removes manually selected ICA components from original (not
% high-pass filtered) data (EEG_ep1).
% curr_resdir: data folder to save new EEG data structure with ICA
% components (EEG_ICA.set) + save selected components (gcompreject.mat)


global ALLCOM ALLEEG CURRENTSET  EEG 

 % data for ICA analysis (it has to be assigned to EEG variable for the eeglab to properly execute ica related functions/GUIs)

%% Reject bad trials
if ~contains(curr_resdir,'LFP')
    
    [EEG_ep1, EEG_ep2] = rej_badtrials(EEG_ep1,EEG_ep2,curr_resdir);
    
end


%% ICA
if ~contains(curr_resdir,'LFP') && length(EEG_ep1.chanlocs)>1
    
    if EEG.trials~=1
        ica_setnm = [curr_resdir filesep 'EEG_ICA.set'];
        crej_nm= [curr_resdir filesep 'gcompreject.mat'];
    else
        ica_setnm = [curr_resdir filesep 'EEG_ICA_continu.set'];
        crej_nm= [curr_resdir filesep 'gcompreject_continu.mat'];
    end
    
    if exist(ica_setnm)==2; ifica = 1; else; ifica = 0; end;
        
    if ifica==0
        
        % get rank of data
        
        % curr_rank = rank(reshape(EEG_ep2.data,[size(EEG_ep2.data,1),size(EEG_ep2.data,2)*size(EEG_ep2.data,3)]));
        
        EEG = EEG_ep2;
        EEG_ICA = pop_runica(EEG, 'icatype', 'fastica');
        pop_saveset(EEG_ICA,ica_setnm);
        
    else
        EEG_ICA = pop_loadset(ica_setnm);
    end
    EEG = EEG_ICA; % it has to be assigned to EEG variable for the eeglab to properly execute ica related functions/GUIs
    
    
    % Label components to reject
    if exist(crej_nm)~=2
        try
            pop_eegplot( EEG, 0, 1, 1,[],'dispchans',20);
            EEG= pop_selectcomps(EEG, 1:20 );
        catch
            fprintf('ICA gone wild.\n')
            close(gcf); close(gcf);
            %continue
        end
        input('Select ICs to reject, if ready, press any key.\n');
    end
    %% Remove selected ICA components from original EEG
    EEG_ep1 = applyica(EEG,EEG_ep1,curr_resdir,crej_nm);
    
end


%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function EEG_ep1 = applyica(EEG,EEG_ep1,curr_resdir,crej_nm)

if exist(crej_nm)~=2
    gcompreject = EEG.reject.gcompreject;
    save(crej_nm,'gcompreject');
else
    load(crej_nm);
end

% Apply ICA for "minimally" filtered data (dataset 1)
EEG_ep1.reject = EEG.reject; EEG_ep1.icawinv = EEG.icawinv;
EEG_ep1.icasphere = EEG.icasphere; EEG_ep1.icaweights = EEG.icaweights; EEG_ep1.icachansind = EEG.icachansind;
EEG_ep1.reject.gcompreject = gcompreject;

EEG_ep1 = pop_subcomp(EEG_ep1,[],1,0);
