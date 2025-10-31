function TE = PD_nosyncTE(Results,interp)
%PD_nosyncTE extracts behavioral information of behavioral analysis
%   PD_nosyncTE(Results,ifinterp) extracts behavioral information from
%   RESULTS struct (data provided by the SSRT task algorithm, see Task_SSRT.m). Event timestamps 
%   are NOT SYNCHRONIZED with EEG data, therefore it is suitable only for 
%   behavioral data analysis). If IFINTERP is true, interpolates values failed to save 
%   by Arduino (due to failed sensor detection), based on timestamps saved by MatLab; if false, discards these values.
%
%   Data is stored in TrialEvents struct and it is saved in the session folder (sess2analyse.folder) under the name:
%               'TrialEvents_nosync_stimoff.mat' |
%               'TrialEvents_nosync_stimon.mat' - in case of postoperative data
%               'TrialEvents_nosync.mat' - in case of intraoperative data
%

% Johanna Petra Szabó, 10.2024
% Lendulet Laboratory of Systems Neuroscience
% Institute of Experimental Medicine, Budapest, Hungary
% szabo.johanna@koki.hun-ren.hu




CUE_ARD = Results.Trial.CUE_ARD(Results.Trial.CUE_ARD>=0)/1000;  % Cue
FB_ARD = Results.Trial.FB_ARD(Results.Trial.FB_ARD>=0)/1000;   % Feedback

if ~any(CUE_ARD~=0);
    fprintf('NO Arduino timestamps.'); TE = [];
    return;
end
validTRinx = 1:length(CUE_ARD);

if interp
    %Interpolate missing ts
    CUE_cog = Results.Trial.flipTime;
    [CUE_ARD_interp, cueintx, CUE_ARD_extrap,cueexpx] = interp_missing_vals(CUE_ARD,CUE_cog);
    CUE_ARD(cueintx) =  CUE_ARD_interp;
    CUE_ARD(cueexpx) =  CUE_ARD_extrap;
    
    
    FB_cog = Results.Trial.fbTime;
    [FB_ARD_interp, cueintx, FB_ARD_extrap,cueexpx] = interp_missing_vals(FB_ARD,FB_cog);
    FB_ARD(cueintx) =  FB_ARD_interp;
    FB_ARD(cueexpx) =  FB_ARD_extrap;
    
    
else
    
    % Discard missing Cue & FB ts
    nocue = find(CUE_ARD==0);
    nofb = find(FB_ARD==0);
    notr = unique([nocue nofb]);
    if ~isempty(notr)
        CUE_ARD(notr) = [];
        FB_ARD(notr) = [];
        validTRinx(notr) = [];
    end
    
end

% CUE
TE.StimulusOn = CUE_ARD;


NumTrials = length(CUE_ARD);   % sync cues
TE.isStopTrial = Results.Trial.isStopTrial(validTRinx);   % assigned as STOP trial; reaction time may come before STOP screen

% Feedback
[TE.Correct, TE.Error] = deal(nan(1,NumTrials));
TE.Feedback = FB_ARD;   % interpolate
% codecheckplot(t,Data,'StartIndex',800,'LineEvent',CueFbTTLOff_MER)
% plot(TE.Feedback,0.1*ones(size(TE.Feedback)),'o')
TE.Outcome = Results.Trial.Outcome(validTRinx);   % 1, correct; 2, error; 3, failed stop; 4, succesful stop
TE.Correct(TE.Outcome==1) = 1;  % 1 for correct press, NaN otherwise
TE.Error(TE.Outcome==2) = 1;  % 1 for incorrect press, NaN otherwies

% Stop screen



TE.StopSignal = nan(1,NumTrials);
[TE.RealStopTrial, TE.FailedStopTrial, TE.SuccesfulStopTrial] = deal(nan(1,NumTrials));

inx = Results.Trial.Situation(validTRinx) == 3 | Results.Trial.Situation(validTRinx) == 4;
TE.RealStopTrial(inx) = 1;   % real stop trials: stop signal presented
inx = Results.Trial.Situation(validTRinx) == 3;
TE.FailedStopTrial(inx) = 1;   % failed to stop
inx = Results.Trial.Situation(validTRinx) == 4;
TE.SuccesfulStopTrial(inx) = 1;   % sucesful stop


STOP_ARD = Results.Trial.STOP_ARD(validTRinx)/1000;
stpi = find(STOP_ARD>0);

stpi2 = find(TE.RealStopTrial==1);
realstopinx = intersect(stpi,stpi2);
SP_ARD =  STOP_ARD(realstopinx);

if interp
    %Interpolate missing ts
    SP_cog = Results.Trial.stopTime;
    [SP_ARD_interp, cueintx, SP_ARD_extrap,cueexpx] = interp_missing_vals(SP_ARD,SP_cog);
    SP_ARD(cueintx) =  SP_ARD_interp;
    SP_ARD(cueexpx) =  SP_ARD_extrap;
end
TE.StopSignal(realstopinx) = SP_ARD;  % stop signla presentation time


TE.PlannedStopDelay = Results.Trial.Delay(validTRinx);  % stop screen delay (Matlab)

% Response, Reaction Time



KEY1_ARD = Results.Trial.RT1_ARD(validTRinx)/ 1000;    % 1st key press time
KEY2_ARD = Results.Trial.RT2_ARD(validTRinx)/ 1000;    % 2nd key press time

TE.Cue1 = Results.Trial.Cue1(validTRinx);  % cue: first number
TE.Cue2 = Results.Trial.Cue2(validTRinx);  % cue: second number
TE.Key1 = Results.Trial.Key1(validTRinx);  % button pressed first
TE.Key2 = Results.Trial.Key2(validTRinx);  % button pressed second
TE.KeyPress1 = KEY1_ARD;   % interpolate 1st key press time
TE.KeyPress2 = KEY2_ARD;   % interpolate 2nd key press time
TE.RT1 = Results.Trial.RT1_ARD(validTRinx)/1000;
TE.RT2 = Results.Trial.RT2_ARD(validTRinx)/1000;
TE.Delay = Results.Trial.Delay(validTRinx)/1000;

% Relative times
TE.TrialStart = TE.StimulusOn - Results.Session.ISI;   % Trial Start: white screen foreperiod
% NOTE: TrialsStart not exact, relies on KbWait.m!!! (can be changed if needed)
TE.StimulusOn = TE.StimulusOn - TE.TrialStart;   % all time stamps relative to TrialStart
TE.Feedback = TE.Feedback - TE.TrialStart;
TE.StopSignal = TE.StopSignal - TE.TrialStart;
TE.KeyPress1 = TE.KeyPress1 - TE.TrialStart;
TE.KeyPress2 = TE.KeyPress2 - TE.TrialStart;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save 'TrialEvents' file




close all



function [CUE_ARD_interp, cueintx, CUE_ARD_extrap,cueexpx] = interp_missing_vals(CUE_ARD_init,CUE_cog)
%INTERP_MISSING_VALS interpolates values failed to saved by Arduino (CUE_ARD_INIT), due to 
%   failed sensor detection, based on MatLab timestamps (CUE_COG).
%
% See also SSRT_TASK, PD_NOSYNC_TE.

nocuei = find(CUE_ARD_init==0); cuei = find(CUE_ARD_init~=0);
if any(nocuei)
    
    cog_cue = CUE_cog(cuei);
    cog_nocue =CUE_cog(nocuei);
    CUE_ARD_interp = interp1(cog_cue,CUE_ARD_init(cuei),cog_nocue);
    
    cog_nocue2 = cog_nocue(isnan(CUE_ARD_interp));
    
    cueintx = nocuei(~isnan(CUE_ARD_interp)); % index of interpolated value in CUE_ARD_init!
    cueexpx = nocuei(isnan(CUE_ARD_interp));
    CUE_ARD_interp = CUE_ARD_interp(~isnan(CUE_ARD_interp));
    
    if ~isempty(cog_nocue2)
        CUE_ARD_extrap=  interp1(cog_cue,CUE_ARD_init(cuei),cog_nocue2,'linear','extrap');
    else
        CUE_ARD_extrap = [];
    end
else
    CUE_ARD_interp = []; CUE_ARD_extrap = [];
    cueintx = [];cueexpx = [];
end

