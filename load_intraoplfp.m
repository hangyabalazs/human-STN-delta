function  EEG =  load_intraoplfp(currsess,patnm,side);

lfpfiles = dir(fullfile(currsess,['*MER*.mat']));

chnr = size(lfpfiles,1);
for mi = 1:chnr
    lfpdat = [];
    lfpdat = load(fullfile(lfpfiles(mi).folder,lfpfiles(mi).name));close(gcf);
    if mi==1
        lfpchans = zeros(chnr,length(lfpdat.Data));
        lfptime = lfpdat.t;
        lfpsr = lfpdat.SampFreq;
    end
    lfpchans(mi,:) = lfpdat.Data;
    if ~isequal(lfpdat.t,lfptime)
        fprintf('Different time vector??? %d - %d\n',mi-1, mi);
    end
end


EEG = pop_importdata('dataformat','array','nbchan',chnr,'data',lfpchans,...
    'setname', [patnm '_' side '_stnlfp' ] ,'srate',lfpsr,...
    'subject',patnm,'pnts',length(lfptime),'xmin',0);
for mii = 1:EEG.nbchan
    EEG.chanlocs(mii).labels = ['Ch' num2str(mii)];
end
EEG.ref = '';
EEG.data = double(EEG.data);