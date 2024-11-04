function [TE_index,evinx] = StimOn_stoppart_evinx(Evinxx,event,evty)
%STIMON_STOPPART_EVINX  Epoch indeces of subevent
% [TE_index,evinx] = STIMON_STOPPART_EVINX(Evinxx,event,evty)
%   -Finds epoch indeces in EVINXX struct of EVENT related epochs, corresponsing to EVTY
%   subevent (Failed or Succesful Stop).
%
% Input parameters:
%   EVINXX  struct storing epoch indeces of event epochs in EEG/ LFP data
%           file/ behavioural data file (see SAVE_EVINXX)
%
%   EVENT   char. array of event label
%
%   EVTY    char. array of subevent label
%
% See also: SAVE_EVINXX


% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu

Tinx1 = Evinxx.(event).(event).TE_index;
Tinx2 = Evinxx.StopSignal.(evty).TE_index;
[TE_index, exx, ~] = intersect(Tinx1,Tinx2);

evinx = Evinxx.(event).(event).epoch_index(exx);