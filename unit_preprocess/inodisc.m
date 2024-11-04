function [TimeStamps, WaveForms] = inodisc(data,ts,sr,thr,varargin)
%INODISC   Unit discrimination.
%   [T W] = INODISC(DATA,TS,SR,THR) performs threshold discrimination of
%   continuous timestamped (TS) unit data (DATA) sampled at the rate SR
%   using the specified threshold (THR). Peak times (T, 'TimeStamps') and
%   spike waveforms (W, 'WaveForms') are saved for each tetrode. A 300 us
%   censored period is applied and the larger spike is kept. Time window
%   for waveform data is set to -300 to 600 us.
%
%   INODISC(DATA,TS,SR,THR,DR) saves the results in the specified directory 
%   (DR). If DR is empty, the data is not saved.
%
%   See also OEDISC and READ_INOMED.

%   Balazs Hangya, Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   24-Nov-2018

% Default arguments
prs = inputParser;
addRequired(prs,'data',@isnumeric)   % raw or filtered data
addRequired(prs,'ts',@isnumeric)   % timestamps
addRequired(prs,'sr',@isnumeric)   % sampling rate
addRequired(prs,'thr',@isnumeric)   % discrimination threshold
addOptional(prs,'resdir',['C:' filesep 'Balazs' filesep '_data' filesep 'Intan' filesep 'Tina' filesep],...
    @(s)ischar(s)|isempty(s))   % results directory
addParameter(prs,'Filtering','enable',@(s)ischar(s)|...
    ismember(s,{'disable','enable'}))   % switch for filtering
parse(prs,data,ts,sr,thr,varargin{:})
g = prs.Results;

% File name
savestr0 = ['save(''' g.resdir filesep 'Ch'];

% Sampling rate
nqf = sr / 2;   % Nyquist freq.
deadtime = 0.0005;   % 300 us dead time
dtp = round(deadtime*sr);   % dead time in data points

% Waveform window
win = [-0.0003 0.0006];   % -300 to 600 us
winp = round(win*sr);   % waveform window in data points

% Threshold discrimintaion
if isequal(g.Filtering,'enable')
    [b,a] = butter(3,[300 6000]/nqf,'bandpass');   % Butterworth filter
    unit = filter(b,a,data);  % filter
elseif isequal(g.Filtering,'disable')
    unit = data;
else
    error('inodisc:InputArg','Unsupported input argument for filtering.')
end
[tvdisc, tpeaks] = disc(unit,thr);   % discriminate (spike times)
tdata = unit;   % data from the current tetrode

% Dead time
dtv = diff(tvdisc);   % ISI
dtpk = diff(tpeaks);   % comparison of neighboring peaks
while any(dtv<dtp)
    censor_inx = dtv < dtp;   % ISI < dead time
    peak_comp = dtpk > 0;
    delete_inx1 = censor_inx & peak_comp;
    delete_inx2 = [0 censor_inx & ~peak_comp];
    tvdisc([find(delete_inx1) find(delete_inx2)]) = [];
    tpeaks([find(delete_inx1) find(delete_inx2)]) = [];
    dtv = diff(tvdisc);   % ISI
    dtpk = diff(tpeaks);   % comparison of neighboring peaks
end
tvdisc(tvdisc<=-winp(1)|tvdisc>=size(data,1)-winp(2)) = [];   % we may lose some spikes near the ends of the file

% Waveform
winx = repmat(tvdisc(:)+winp(1),1,sum(abs(winp))) + repmat(0:sum(abs(winp))-1,length(tvdisc),1);
wv = nan(size(winx,1),size(winx,2));   % waveforms: spikes x time
wv(:,:) = tdata(winx);   % waveform data
WaveForms = wv;

% Spike times
spike_times = ts(tvdisc);
TimeStamps = spike_times;
savestr = [savestr0 ''',''WaveForms'',''TimeStamps'');'];

% Save
if ~isempty(g.resdir)
    eval(savestr)
end