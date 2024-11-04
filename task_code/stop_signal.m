function [early, respt, early2, respt2, situation, ARD_time, ARD_time2,...
    STOP_ARD, STOP_MTL, stopTime, MTL_time, MTL_time2] = stop_signal(duration, win, stopdelay, S_time, red, serial, CUE_ARD)
%STOP_SIGNAL
%   [EARLY, RESPT, EARLY2, RESPT2, SITUATION] = STOP_SIGNAL(DURATION, WIN, STOPDELAY, S_TIME, RED)
%   Waits for user keypress for DURATION seconds. If a key specified in the
%   code is pressed, it returns outcome (EARLY) and time of button press
%   (RESPT) relative to the start of the cue presentation (Screen, 'Flip').
%   WIN - Cue screen with the presetted parameters.
%   Returns SITUATION number which will be used for feedback.
%   STOPDELAY - waiting period after stop signal where the patient should
%   withold keypress.
%   S_TIME - Start time of the cue presentation ('Flip')
%   RED - Background color for Stop signal
%
%   Implemented keys:
%       '1, 2, 3' - response keys
%       'Esc' - exit task
%
%   Output variables:
%       EARLY - 1 if 'a number' was pressed; NaN if 'Esc' pressed; 0 otherwise
%       RESPT - time of buttonpress if 'a number' was pressed;
%       0 if keypress was withold, NaN in the case of 'Esc' was pressed
%
%   See also DEFAULT_CYCLE3, LISTENKEYS, GETRESPONSE and SUSATTN_MAIN.
%
%   Edit log: BH 7/25/12, TL 4/28/16, TL 6/23/16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define keys
Start = S_time;   % set time 0 (for reaction time)
global Sinx
Sinx.case = 0;   % Stop cycle index 1 - if Stop cycle, 0 - if Default cycle
Sinx.buttoncode = [];
Sinx.ARD_time = [];
Sinx.MTL_time = [];
global RectSize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Listen to key press
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
buttoncode = [];
disp('look for early keypress')
if duration == stopdelay
    duration = 1.99;
end
[buttoncode, ARD_time, MTL_time] = ser_read(serial, duration, Start, win, 2);   % Look for early keypress before stop presentation


if ~isempty(buttoncode)   % EARLY KEYPRESS BEFORE STOP SIGNAL, continue as default cycle
    disp('EARLY KEYPRESS1')
    early = buttoncode;   % Store the value of the pressed key
    respt = ARD_time - CUE_ARD;   % RESPONSE TIME: Elapsed time between Screen 'Flip' and keypress
    buttoncode = [];
    [buttoncode, ARD_time2, MTL_time2] = ser_read(serial, Inf, Start, win, 2);
    disp('EARLY KEYPRESS2')
    early2 = buttoncode;   % Store the value of the second keypress
    respt2 = ARD_time2 - CUE_ARD;   % RESPONSE TIME: Elapsed time between Screen 'Flip' and second keypress
    situation = 2;   % Keypress occured before stop signal, no stop signal presentation, run as a default cycle
    STOP_ARD = 0;
    STOP_MTL = 0;
    stopTime = 0;
    
else   % NO KEYPRESS BEFORE STOP SIGNAL
    disp('No keypress before stop')
    Screen('FillRect',win, red);   % set the background color for STOP signal
    DrawFormattedText(win, 'STOP', 'center', 'center', 250);   % STOP signal
    Screen('Fillrect', win, [0 0 0], [0 0 RectSize]);   % Black square on the top left corner for the photosensor
    
    % STOP SIGNAL PRESENTATION
    [STOP_ARD, STOP_MTL, stopTime] = flipSens(serial, win, 4);
    
    if Sinx.case == 1   % Keypress right before Stop presentation
        disp('SINX1: ALMOST FLIPPED')
        buttoncode = Sinx.buttoncode;
        ARD_time = Sinx.ARD_time;
        MTL_time = Sinx.MTL_time;
        disp('EARLY KEYPRESS1')
        early = buttoncode;   % Store the value of the pressed key
        respt = ARD_time - CUE_ARD;   % RESPONSE TIME: Elapsed time between Screen 'Flip' and keypress
        buttoncode = [];
        [buttoncode, ARD_time2, MTL_time2] = ser_read(serial, Inf, Start, win, 2);
        disp('EARLY KEYPRESS2')
        early2 = buttoncode;   % Store the value of the second keypress
        respt2 = ARD_time2 - CUE_ARD;   % RESPONSE TIME: Elapsed time between Screen 'Flip' and second keypress
        situation = 2;   % Keypress occured before stop signal, no stop signal presentation, run as a default cycle
        STOP_ARD = 0;
        STOP_MTL = 0;
        stopTime = 0;
        Sinx.case = 0;   % Stop cycle index 1 - if Stop cycle, 0 - if Default cycle
        Sinx.buttoncode = [];
        Sinx.ARD_time = [];
        Sinx.MTL_time = [];
    elseif Sinx.case == 2   % Keypress after Stop, before sensor detection
        buttoncode = Sinx.buttoncode(1);
        disp(['SINX2: Buttonpress after flip, before sensor detection: ' num2str(buttoncode)])
        ARD_time = Sinx.ARD_time(1);
        MTL_time = Sinx.MTL_time(1);
        early = buttoncode;
        respt = ARD_time - CUE_ARD;
        if isnan(stopTime)   % STOP NOT DETECTED, CONTINUE WITHOUT TIMESTAMP
            disp('!!! NO STOP FLIP TIMESTAMP !!!')
            if length(Sinx.buttoncode) > 1
                buttoncode = Sinx.buttoncode(2);
                ARD_time2 = Sinx.ARD_time(2);
                MTL_time2 = Sinx.MTL_time(2);
            else
                buttoncode = Sinx.buttoncode(1);
                ARD_time2 = Sinx.ARD_time(1);
                MTL_time = Sinx.MTL_time(1);
            end
            early2 = buttoncode;
            respt2 = ARD_time - CUE_ARD;
            stopTime = Sinx.flipTime;   % ONLY MTL FLIPTIME IS STORED FOR ESTIMATION
            STOP_ARD = 0;
            STOP_MTL = NaN;   % SET IT TO NAN FOR INDICATING THE TRIAL WHEN IT HAPPENED
        else   %  KEYPRESS AFTER FLIP BEFORE SENSOR DETECTION, LOOK FOR 2ND PRESS
            [buttoncode, ARD_time2, MTL_time2] = ser_read(serial, stopdelay, stopTime, win, 2);   % Looking for 2nd keypress
            if ~isempty(buttoncode)   % if second keypress occured saves it
                disp('SINX2: Also second keypress occured after stop')
                early2 = buttoncode;   % Store the value of the second keypress
                respt2 = ARD_time2 - CUE_ARD;   % RESPONSE TIME: Elapsed time between Screen 'Flip' and second keypress
            else   % If only one keypress occured after STOP signal, it sets the second keypress value and RT to 0
                early2 = 0;
                respt2 = 0;
                ARD_time2 = 0;
                MTL_time2 = 0;
                buttoncode = early;
            end
            Sinx.case = 0;   % Stop cycle index 1 - if Stop cycle, 0 - if Default cycle
            Sinx.buttoncode = [];
            Sinx.ARD_time = [];
            Sinx.MTL_time = [];
            Sinx.flipTime = [];
        end
        situation = 3;   % Participant could not withold keypress after stop signal
        
    elseif Sinx.case == 0
        buttoncode = [];
        disp('No keypress during flip, look for press after stop presentation')
        [buttoncode, ARD_time, MTL_time] = ser_read(serial, stopdelay, stopTime, win, 5);   % Looking for keypress for 2s
        if ~isempty(buttoncode)   % KEYPRESS AFTER STOP SIGNAL
            early = buttoncode;
            respt = ARD_time - CUE_ARD;
            situation = 3;   % Participant could not withold keypress after stop signal
            buttoncode = [];
            disp('Keypress after stop')
            [buttoncode, ARD_time2, MTL_time2] = ser_read(serial, stopdelay, stopTime, win, 2);   % Looking for 2nd keypress
            if ~isempty(buttoncode)   % if second keypress occured saves it
                disp('Also second keypress occured after stop')
                early2 = buttoncode;   % Store the value of the second keypress
                respt2 = ARD_time2 - CUE_ARD;   % RESPONSE TIME: Elapsed time between Screen 'Flip' and second keypress
            else   % If only one keypress occured after STOP signal, it sets the second keypress value and RT to 0
                early2 = 0;
                respt2 = 0;
                ARD_time2 = 0;
                MTL_time2 = 0;
                buttoncode = early;
            end
        end
    end
    
    if isempty(buttoncode)   % NO KEYPRESS AFTER STOP SIGNAL, all parameters set to 0
        early = NaN;
        respt = NaN;
        early2 = NaN;
        respt2 = NaN;
        ARD_time = NaN;
        ARD_time2 = NaN;
        MTL_time = NaN;
        MTL_time2 = NaN;
        situation = 4;   % Succesfull keypress inhibition
        disp('Succesfull keypress inhibition')
    end
end

if buttoncode == 4   % ESC PRESSED, EXITING THE LOOP, sets all parameters to NaN's
    early = NaN;
    respt = NaN;
    early2 = NaN;
    respt2 = NaN;
    ARD_time = NaN;
    ARD_time2 = NaN;
    MTL_time = NaN;
    MTL_time2 = NaN;
    situation = 5;   % ESC pressed, ending the trial
    disp('ESC pressed')
end
