function [early, respt, early2, respt2, situation, ARD_time,...
    ARD_time2, MTL_time, MTL_time2] = default_cycle(flipTime, ser, win)
%DEFAULT CYCLE
%   [EARLY, RESPT, EARLY2, RESPT2, SITUATION, ARD_time, ARD_time2]...
%   = DEFAULT_CYCLE(FLIPTIME, SER, WIN)
%   Waits for user keypresses. If a keypress occures, it returns it's code
%   (EARLY) and time of button press (RESPT) relative to the start of the
%   cue presentation (Screen, 'Flip'), detected by the sensor.
%   Returns SITUATION number which will be used for feedback.
%   FLIPTIME - Start time of the cue presentation ('Flip')
%
%   Implemented keys:
%       '1, 2, 3' - response keys
%       'terminating button (Esc)' - exit task
%
%
%   See also STOP_SIGNAL, FLIPSENS, SER_READ and FGETL_NONBLOCKING.
%
%   Edit log: BH 7/25/12, TL 4/28/16, TL 6/23/16, TL 25/01/17

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Listen to key press
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Start = flipTime;   % Flip execution time
global Sinx
Sinx.case = 0;   % Stop cycle index 1 - if Stop cycle, 0 - if Default cycle
Sinx.buttoncode = [];
Sinx.ARD_time = [];
Sinx.MTL_time = [];
buttoncode = [];
% Checking the buttons till keypress occures
disp('DEFAULT KEY1')
[buttoncode, ARD_time, MTL_time] = ser_read(ser, Inf, Start, win, 2);


if buttoncode == 4 % ESC pressed, set all parameters to NaN
    early = NaN;
    respt = NaN;
    early2 = NaN;
    respt2 = NaN;
    ARD_time = NaN;
    ARD_time2 = NaN;
    MTL_time = NaN;
    MTL_time2 = NaN;
    situation = 5;   % ESC pressed, exit the loop, end of trial
else
    early = buttoncode;
    respt = ARD_time - Start;   % RT: Time between 'Flip' and keypress
    buttoncode = [];
end

% Checking for second keypress
disp('DEFAULT KEY2')
if isempty(buttoncode)
    [buttoncode, ARD_time2, MTL_time2] = ser_read(ser, Inf, Start, win, 2);
end


if buttoncode == 4 % ESC pressed, set all parameters to NaN
    early = NaN;
    respt = NaN;
    early2 = NaN;
    respt2 = NaN;
    ARD_time = NaN;
    ARD_time2 = NaN;
    MTL_time = NaN;
    MTL_time2 = NaN;
    situation = 5;   % ESC pressed, exit the loop, end of trial
else
    % Default situation, when 2 keypress occured after the cue presentation
    early2 = buttoncode;
    respt2 = ARD_time2 - Start;   % RT: Time between 'Flip' and keypress
    buttoncode = [];
    situation = 1;
end


