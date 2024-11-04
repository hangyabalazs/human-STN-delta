function saveDATasMER(datadir)
%SAVEDATASMER converts inomed '*.dat' file to '*.mat' file
% SAVEDATASMER(datadir) reads the manually chosen '*.dat' file from directory DATADIR, plots the
% recorded signal and converts it to a '*.mat' file.

% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% hnn

cd(datadir)
[mer_file, mer_filedir] = uigetfile('*EMG*.dat', ...
    'Choose a DAT file!', 'MultiSelect', 'on');

if iscell(mer_file)
    NumChannels = length(mer_file);% number of channels
else
    NumChannels =1;
end

for iT = 1:NumChannels
    
    if NumChannels>1
        fileName = mer_file{iT};
    else
        fileName = mer_file;
    end
    
    fid = fopen(fileName);
    
    MF.Header=fread(fid,80,'*char')';
    MF.SampFreq=double(fread(fid,1,'*int')');             % 80..81
    
    MF.SiteNr=double(fread(fid,1,'*int')');
    MF.KanalNr=double(fread(fid,1,'*int')');
    MF.op_id=double(fread(fid,1,'*int')');
    MF.SiteID=double(fread(fid,1,'*int')');
    MF.MaxYValue=double(fread(fid,1,'*int')');
    MF.typ=double(fread(fid,1,'*int')');
    MF.Multiplexing=double(fread(fid,1,'*int')');
    MF.EMGChannels=double(fread(fid,1,'*int')');
    
    
    MF.BitsPerValue=double(fread(fid,1,'*int')');
    MF.fcHPInHz=double(fread(fid,1,'*int')');
    MF.fcLPInHz=double(fread(fid,1,'*int')');
    MF.vu=double(fread(fid,1,'*int')');
    MF.fcHPInHzSW=double(fread(fid,1,'*int')');
    MF.fcLPInHzSW=double(fread(fid,1,'*int')');
    MF.bitResolution=double(fread(fid,1,'*double')');
    MF.versionNrMER=char(fread(fid,32,'*char')');
    
    MF.Data=(double(fread(fid,20000000,'*ushort'))-(MF.MaxYValue/2))*(0.4/MF.MaxYValue);
    MF.Data = (MF.Data-mean(MF.Data))*10^4;
    MF.t=0:(1/MF.SampFreq):((length(MF.Data)-1)/MF.SampFreq);
    graph=figure;
    plot(MF.t,MF.Data);
    title('Data');
    xlabel('Time in s');
    ylabel('Voltage in mV');
    
    
    fileName=fileName(1:length(fileName)-4);
    saveas(graph,[fileName '_plot.fig']);
    saveas(graph,[fileName '_plot'], 'jpg');
    
    fclose(fid);
    
    
    save([fileName '.mat'], '-struct','MF');
    
end