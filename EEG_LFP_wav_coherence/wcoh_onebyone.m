function wcoh_onebyone(sess2analyse,event,subevent,isfig,freqs,downsamp,alpha,plot_win,baseline_win,chanmean)
% WCOH_ONEBYONE  Mean-squared wavelet coherence (MSWC) map between intraop EEG and LFP.
%   EEG_LFP_WCOH_PD calculates the MSWC map aligned to EVENT within FREQS frequency range 
%   patient-by-patient, including trials partitioned according to SUBEVENT.
%   If DOWNSAMP, number of partitioned trials is downsampled to have equal
%   nr. of trials in each partitioned set. 
%   Significant changes relative to BASELINE_WIN time window are tested
%   with permutation test with ALPHA level, using FDR correction.
%   If ISFIG, draw MSWC map for each patient within PLOT_WIN time window (sec rel. to event start).
%   If CHANMEAN, LFP is calculated as the average of recording channels.
% 
% Input parameters:
%     SESS2ANALYSE      struct containing all necessary information (name of patient, side
%                       of experiment, tag of condition, session folder path) of
%                       session data that need to be analysed (see getdata2analyse)
%
% See also: WCOHERENCE, STATCONDFIELDTRIP, GET_PATIENT_WCOH
%
% Johanna Petra Szabó, 05.2025
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

global figdir_pd

mcorrect = 'cluster';
if isempty(subevent)
    figdir = fullfile(figdir_pd,'introp_EEG_LFP','WCoh',event);
    subevent= event;
else
    figdir = fullfile(figdir_pd,'introp_EEG_LFP','WCoh',[event '_' subevent]);
end
if ~isfolder(figdir); mkdir(figdir); end;

time = -2 :1/250: 2; time = time(1:end-1);
plotlim_inx = dsearchn(time',plot_win'); plotinx = plotlim_inx(1):plotlim_inx(2);
pltime = time(plotinx);
if ~isempty(baseline_win)
    baslim_inx = dsearchn(time',baseline_win'); basinx = baslim_inx(1):baslim_inx(2);
else
    basinx = plotinx;
end

wcoh_structname = fullfile(figdir, ['WavCH_WIN' num2str(plot_win) '_BAS' num2str(baseline_win)]);
stat_structname = fullfile(figdir, ['WavCH_STAT_WIN' num2str(plot_win) '_BAS2' num2str(baseline_win)]);

if exist([wcoh_structname '.mat'])~=2
    WavCH  = struct;
else
    load([wcoh_structname '.mat']);
end

if exist([stat_structname '.mat'])~=2
    WavCH_STAT = struct;
else
    load([stat_structname '.mat']);
end

eeg_chan = 'F4';
for k = 1:length(sess2analyse)
    
    patnm = sess2analyse(k).patient;
    side = sess2analyse(k).side;
    fprintf('%s %s ...\n',patnm,side);
    
    %% EEG
    
    % Load EEG data
    curr_resdir_eeg = sess2analyse(k).folder;
    
    
    [wcoh,phlag,f,coi,chnames,~] = get_patient_wcoh(patnm,side,curr_resdir_eeg,event,subevent,downsamp);
    
    if isempty(wcoh);
        continue;
    end
    
    if chanmean
        wcoh000 = cell(1,length(chnames));
        for ci = 1:length(chnames);
            wcoh000{ci} = cat(3,wcoh{:,ci});
        end
        wcoh_chan = nanmean( cat(4,wcoh000{:}) ,4);
        chnames = {'chanmean'};
    end
    
    ff = f{1}; coicoi = coi{1};
    
    for ci = 1:length(chnames);
        lfp_chan = chnames{ci};
        
        % Average wavelet coherence over epochs
        if ~chanmean
            wcoh_chan = cat(3,wcoh{:,ci});
        end
        
        if ~isempty(baseline_win)
            wcoh_bas_avg = repmat( nanmean(wcoh_chan(:,basinx,:),[2 3]) ,[1 size(wcoh_chan,2) 1]);
            wcoh_bas_std = repmat( std(wcoh_chan(:,basinx,:),[],[2 3],'omitnan') ,[1 size(wcoh_chan,2) 1]);
            wcoh_chanavg = (nanmean(wcoh_chan,3)-wcoh_bas_avg)./wcoh_bas_std;
            
%             wcoh_bas_avgx = repmat( nanmean(wcoh_chan(:,basinx,:),2) ,[1 size(wcoh_chan,2) 1]);
%             wcoh_bas_stdx = repmat( std(wcoh_chan(:,basinx,:),[],2,'omitnan') ,[1 size(wcoh_chan,2) 1]);
%             wcoh_chanN = (wcoh_chan - wcoh_bas_avgx)./ wcoh_bas_stdx;
            
        else
            wcoh_chanavg = nanmean(wcoh_chan,3);
            wcoh_chanN = wcoh_chan;
        end
        
%             wcoh_chanavg = nanmean(wcoh_chan,3);
        % Average phase lag
        phlag_chan = cat(3,phlag{:,ci});
        phlag_chanavg = angle(  nanmean( exp(1i*phlag_chan) ,3 )  );
%         ITphlag_chanavg = abs(  nanmean( exp(1i*phlag_chan) ,3 )  );
               
        % Stat
        
          
        f_ind = find(ff>freqs(1)&ff<freqs(2));
        fff = ff(f_ind);
        
        
        
        if ~isnan(alpha)
%              [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(wcoh_chan(:,plotinx,:),ff,alpha,2000,false,'fdr');
            [exactp_ersp,maskersp,alphafdr] = boostat_eeglab_J(wcoh_chan(:,plotinx,:),ff,alpha,1000,false,'fdr',[],1:length(basinx));
                       
            
            
            % Save  coherence map STATs
            WavCH_STAT.(patnm).(side).([eeg_chan '_' lfp_chan]).mask_ersp = maskersp;
            WavCH_STAT.(patnm).(side).([eeg_chan '_' lfp_chan]).p_ersp = exactp_ersp;
%             WavCH_STAT.(patnm).(side).([eeg_chan '_' lfp_chan]).diffmap = diffmap;
            WavCH_STAT.(patnm).(side).([eeg_chan '_' lfp_chan]).alphafdr = alphafdr;
            WavCH_STAT.(patnm).(side).([eeg_chan '_' lfp_chan]).alpha = alpha;
            
            
            if ~downsamp
                save([stat_structname '.mat'],'WavCH_STAT');
            else
                save([stat_structname '_downsamp.mat'],'WavCH_STAT');
            end
            
            
        end
      
        
        if isfig
            % TIME-FREQ. W-COH PLOT
            fig = figure;
            wcoh_fig(wcoh_chanavg(f_ind,plotinx),phlag_chanavg(f_ind,plotinx),pltime,fff,coicoi(plotinx),[-.4 .4]);
            if ~isnan(alpha);
                hold on;
                
                zmap = norminv(exactp_ersp);
                h= bootstatFDR_clustercorr(zmap(f_ind,:),maskersp(f_ind,:),pow2(log2(fff)),'intraop','EEG_LFP','Fp1',pltime,log2(fff))
                %                 contour(pltime,log2(fff),maskersp,'Color','red');
            end
            arrayfun(@(x) set(x,'Color','r'), h(~isnan(h)),'UniformOutput',0);
            % Save plot
            fnm= [patnm '_' side '_' event '_' eeg_chan '_' lfp_chan '_F' num2str(freqs) '_WIN' num2str(plot_win) '_BAS2' num2str(baseline_win)];
%             fnm= [patnm '_' event '_' eeg_chan '_' lfp_chan '_F' num2str(freqs) '_WIN' num2str(plot_win) '_BAS' num2str(baseline_win)];
            saveas(gcf,fullfile(figdir,[fnm '.jpg']))
            saveas(gcf,fullfile(figdir,[fnm '.fig']))
%             saveas(gcf,fullfile(figdir,[fnm '.pdf']))
            close(gcf);
            
%             % PHASE LAG consistency (resultant vector length)
%             fig = figure;
%             wcoh_fig(ITphlag_chanavg(f_ind,:),phlag_chanavg(f_ind,:),time,fff,coicoi);
%             colormap(pink)
%             caxis([0 0.3])
%             
%             % Save plot
%             fnm= [patnm '_' event '_' eeg_chan '_' lfp_chan '_F' num2str(freqs) '_meanPhLag'];
%             saveas(gcf,fullfile(figdir,[fnm '.jpg']))
%             saveas(gcf,fullfile(figdir,[fnm '.fig']))
%             close(gcf);
        end
        
        
        
        % Save averaged coherence map
        WavCH.(patnm).(side).([eeg_chan '_' lfp_chan]) = wcoh_chanavg;
        
        if ~downsamp
            save([wcoh_structname '.mat'],'WavCH');
        else
            save([wcoh_structname '_downsamp.mat'],'WavCH');
        end
        save(fullfile(figdir,'WavCH_f.mat'),'ff');
        save(fullfile(figdir,'WavCH_coi.mat'),'coicoi');
    end
end
end