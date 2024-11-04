function Task_SSRT
%TASK_SSRT   Stop signal reaction time task for Parkinson patients

%%%%% AIMS %%%%%
% The task uses a buttonbox which has 3 buttons, a monitor, a
% lightsensor, and a trigger output. The participant has to press two
% buttons in the same order as it is presented in the cue. After the
% keypress there is a feedback screen giving positive or negative
% reinforcement. In a few cases after the cue a stop signal is presented.
% Participants are tested if they can withold keypressing based on the stop
% latency. The reaction times are constantly monitored and stop latencies
% are adjusted trial-by-trial.
% Button presses and cue presentation are detected and timestamped by
% an Arduino. Arduino code, and serial read done by Andras Szell.

% 02/02/2016 created by Tamas Laszlovszky

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

% Initialize pseudo random stream
rng('shuffle');
dbstop if error;

% Select folder
Screen('Preference', 'SkipSyncTests', 2);  % Sync problem with Windows 10, force to run Psychtoolbox
results_folder = 'C:\Users\ntadmin\Documents\MATLAB\_PD_DATA';
cd(results_folder);   % Set current folder
% Identify the current session type
fprintf('Adja meg a meres tipusat: \n1 - Preop \n2 - Intraop \n3 - Postop \n')
currType = str2num(input('Jelenlegi meres kodja?\n','s'));

recordTypes = {'preop', 'intraop', 'postop'};
recordType = recordTypes{currType};

% Task parameters
switch recordType
    case 'preop'
        StimTypes = {'1 2' '1 3' '3 2' '2 1'};   % Number combinations for cues (misses 2-3 and 3-1
    case 'intraop'
        StimTypes = {'1 2' '1 3' '2 3' '3 1' '3 2' '2 1'};   % Number combinations for cues
    case 'postop'
        StimTypes = {'1 2' '1 3' '2 3' '3 1' '3 2' '2 1'};   % Number combinations for cues
end

% Select side
fprintf('Melyik oldalt tesztelik?: \n1 - Bal \n2 - Jobb \n')
currSide = str2num(input('Oldal kodja?\n','s'));
sides = {'left', 'right'};
side = sides{currSide};


% Identify the Patient
trialNum = input('Beteg bevonasi szama?\n','s');
if isempty(trialNum)
    error('Nem adott meg beteg azonosítót. A program leállt.')
end
letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';

% Create unique filename for results
next = 1;
results_fnm = results_folder;
while exist(results_fnm,'file')
    results_fnm = [results_folder filesep 'n' trialNum '_' recordType '_' side '_' datestr(now,'yymmdd') letters(next) '.mat'];
    next = next + 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTING THE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ISI = 0.5;  % foreperiod (seconds)
background = [255 255 255];   % Default background color
green = [0 214 0];   % Hit background color
red = [214 0 0];   % False background color
stopdelay = 2;   % 2 second waiting period for stop signal inhibition
[width, height]=Screen('DisplaySize', 2);
CompScreen = get(0,'MonitorPositions');% Find out the size of this computer screen
CompScreen = CompScreen(2,:);
global RectSize
% RectSize = [CompScreen(3)*0.05 CompScreen(3)*0.05];
%RectSize = [CompScreen(1)*-0.04 (CompScreen(4)-CompScreen(2))*0.12];
RectSize = [70, 40];
minlatency = 0.1;   % min stop signal latency is 0.1 s
maxlatency = 2;   % max stop signal latency
NumTrials = 1000;   % maximum number of trials (for preallocation)
StimNums1 = NaN(1, NumTrials);   % Cues first number
StimNums2 = NaN(1, NumTrials);   % Cues second number
StimList = StimTypes(randi([1 length(StimTypes)],1,NumTrials));  % randomize trial order
for i = 1:NumTrials   %store them as numbers too
    StimNums1(i) = str2num(StimList{i}(1));
    StimNums2(i) = str2num(StimList{i}(end));
end
isStopTrial = rand(1,NumTrials) < 0.6;   % draw stop trials (40%)
isStopTrial(1:5) = false;   % no stop signal in the first 5 trials
for i = 1:NumTrials-1
   if isStopTrial(i) == 1
       if isStopTrial(i+1) == 1
           isStopTrial(i+1) = 0; % no consecutive Stop trials
       end
   end
end
global backup   % global variable for storing every bit coming from the serial buffer
backup = [];
global backup_buff   % temporary backup in case of buffer overfill
backup_buff = [];
global Sinx   % Stop Index: Pre-correction variable for sensor and button data
Sinx.buttoncode = [];
Sinx.ARD_time = [];
Sinx.MTL_time = [];
Sinx.flipTime = [];
Sinx.case = [];
global cue; % Storing current cue combination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
Results.Session.ISI = ISI;
Results.Session.NumTrials = NumTrials;
Results.Session.Cycle = NaN;
Results.Session.trialType = recordType;
Results.Session.Side = side;
Results.Session.Subject = trialNum;
Results.Session.MinLatency = minlatency;
Results.Session.MaxLatency = maxlatency;
Results.Session.ScreenSize = CompScreen;
Results.Session.StopDelay = stopdelay;

% Preallocate
Results.Trial.Key1 = NaN(1, NumTrials);
Results.Trial.Key2 = NaN(1, NumTrials);
Results.Trial.RT1 = zeros(1, NumTrials);
Results.Trial.RT2 = zeros(1, NumTrials);
Results.Trial.RT1_MTL = zeros(1, NumTrials);
Results.Trial.RT2_MTL = zeros(1, NumTrials);
Results.Trial.MTL_time = zeros(1, NumTrials);
Results.Trial.MTL_time2 = zeros(1, NumTrials);
Results.Trial.Outcome = NaN(1, NumTrials);
Results.Trial.Cue1 = StimNums1;
Results.Trial.Cue2 = StimNums2;
Results.Trial.isStopTrial = isStopTrial;
Results.Trial.Delay = NaN(1, NumTrials);
Results.Trial.Situation = NaN(1, NumTrials);
Results.Trial.TrialLength = [];
Results.Trial.outcome = NaN(1, NumTrials );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Present Instruction Screen
% win = Screen('OpenWindow',2, [900 900 1000], CompScreen); % Full sized screen
win = Screen('OpenWindow',2, [900 900 1000], [-1365,1, 1,768]);% Smallscreen used for testing
      white = WhiteIndex(win);   % white for current screen depth
Results.Session.WhiteI = white;

% LightSensor configuration
Screen('FillRect', win, white);   % Create white screen
Screen('Fillrect', win, [0 0 0], [0 0 RectSize]);   % Black square on the top left corner for the photosensor
Screen(win,'Flip'); % present to the screen. This is the command to actually present whatever you have made 'win'
KbWait;
WaitSecs(ISI);
Screen('FillRect',win, white);   % fill the screen with white
Screen(win,'Flip'); % present to the screen
KbWait;
WaitSecs(ISI);

% Instruction screen
Screen('TextSize',win, 60);
Screen('TextFont',win,'Courier New');
Screen('TextStyle', win, 1);
[~, screenYpixels] = Screen('WindowSize', win);   % Compute the screen center
DrawFormattedText(win, 'Parkinson SSRT Teszt', 'center', screenYpixels * 0.1, 0);
DrawFormattedText(win, 'Üdvözöljük!', 'center', screenYpixels * 0.4, 0);
DrawFormattedText(win, 'Amint készen áll, indítjuk a tesztet.', 'center', screenYpixels * 0.6, 0);
Screen('TextSize',win, 60);
Screen(win,'Flip');
TrialLength = GetSecs;   % Start measuring overall length of the experiment
WaitSecs(ISI); % this delay avoids participants accidently pressing too quickly
KbWait;
WaitSecs(ISI);


% Initialize serial comport
delete(instrfindall())
ser = serial('COM3', 'Baudrate', 115200);
fopen(ser);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN TRIALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Trial loop
for iT = 1:NumTrials    % runs until maximum number of trials is reached or pressing ESC
    disp(['Trial Number: ' num2str(iT)])
    % Present cue
    cue = StimList{iT};   % current stimulus
    Screen('TextSize',win, 300);
    DrawFormattedText(win, cue, 'center', 'center', 0);          % Cue presentation
    
    % SELECTING DEFAULT OR STOP CYCLE
     if isStopTrial(iT)   % Enters only after the first 5 cycles and in case of isStopTrial = 1
         % STOP LOOP
        avgRTp7 = mean(Results.Trial.RT1(iT-5:iT-1)) * maxlatency;   % Calculates avg reaction time in the last 5 cycles * maxlatency
        delay = (avgRTp7 - minlatency).*rand + minlatency;   % stop signal latency in a range based on min and max latency
        if isnan(delay) || delay>100 % If no flip detection switch to MTL_time based reaction time latency
            disp('SWITCH TO MTL TIMESTAMPS FOR DELAY')
            avgRTp7_MTL = nanmean(Results.Trial.RT1_MTL(iT-5:iT-1)) * maxlatency;   % Calculates avg reaction time * maxlatency
            delay_MTL = (avgRTp7_MTL - minlatency).*rand + minlatency;   % stop signal latency in a range based on min and max latency
            delay = delay_MTL;
            disp(num2str(delay))
        end
        Screen('Fillrect', win, [0 0 0], [0 0 RectSize]);   % Black square on the top left corner for the photosensor
        [CUE_ARD, CUE_MTL, flipTime] = flipSens(ser, win, 1);   % Cue presentation, sensor detection
        [early, respt, early2, respt2, situation, ARD_time, ARD_time2,...
            STOP_ARD, STOP_MTL, stopTime, MTL_time, MTL_time2] = stop_signal(delay, win,...
            stopdelay, flipTime, red, ser, CUE_ARD);   % ENTERING STOP CYCLE
     else
        % DEFAULT LOOP
        situation = 1;   % Default cycle when no stop signal is presented
        delay = 0;   % Reset delay
        Screen('Fillrect', win, [0 0 0], [0 0 RectSize]);   % Black square on the top left corner for the photosensor
        [CUE_ARD, CUE_MTL, flipTime] = flipSens(ser,win, 1);   % Cue presentation, sensor detection
        [early, respt, early2, respt2, situation, ARD_time, ARD_time2, MTL_time, MTL_time2] = default_cycle(CUE_ARD, ser, win);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUATING THE RESPONSES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    switch situation
        case 1   % Default cycle
            if early == StimNums1(iT) && early2 == StimNums2(iT)
                hit = 1;
                feedback = 'Jó válasz';
                background = green;
                outcome = 1;
            else
                hit = 2;
                feedback = 'Rossz válasz';
                background = red;
                outcome = 2;
            end
        case 2   % Early keypress before stop signal, continues as a default cycle
            if early == StimNums1(iT) && early2 == StimNums2(iT)
                hit = 1;
                feedback = 'Jó válasz';
                background = green;
                outcome = 3;
            else
                hit = 2;
                feedback = 'Rossz válasz';
                background = red;
                outcome = 4;
            end
        case 3   % Keypress occured after Stop signal presentation
            feedback = 'Téves gombnyomás \nSTOP jelzés után';
            hit = 3;
            background = red;
            outcome = 5;
        case 4   % Successfull stop signal inhibition
            feedback = 'Sikeres gombnyomás \nmegállítás';
            hit = 4;
            background = green;
            outcome = 6;
            
        case 5   % ESC pressed, End of Trial
            break
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FEEDBACK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Feedback
    Screen('FillRect',win, white);   % fill the screen with white
    Screen('Flip', win);  % present to the screen
    WaitSecs(ISI);
    Screen('FillRect',win, background);   % set background color according to feedback
    Screen('TextSize',win, 90);
    DrawFormattedText(win, feedback, 'center', 'center', 250);          % Feedback presentation
    Screen('Fillrect', win, [0 0 0], [0 0 RectSize]);   % Black square on the top left corner for the photosensor
    [FB_ARD, FB_MTL, fbTime] = flipSens(ser, win, 3);   % Feedback presentation, sensor detection
    WaitSecs(2);
    Screen('FillRect',win, white);   % fill the screen with white
    Screen('Flip', win);   % present to the screen
    WaitSecs(ISI);   % next trial foreperiod
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STORING THE DATA OF THE ACTUAL CYCLE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Results.Trial.Key1(iT) = early;
    Results.Trial.Key2(iT) = early2;
    Results.Trial.RT1(iT) = respt/1000;
    Results.Trial.RT2(iT) = respt2/1000;
    if isempty(MTL_time)
        MTL_time = flipTime;
    end
    if isempty(MTL_time2)
        MTL_time2 = NaN;
    end
    Results.Trial.RT1_MTL(iT) = MTL_time - flipTime;
    Results.Trial.RT2_MTL(iT) = MTL_time2 - flipTime;
    Results.Trial.RT1_ARD(iT) = ARD_time;
    Results.Trial.RT2_ARD(iT) = ARD_time2;
    Results.Trial.MTL_time(iT) = MTL_time;
    Results.Trial.MTL_time2(iT) = MTL_time2;
    Results.Trial.Outcome(iT) = hit;
    Results.Trial.Delay(iT) = delay;
    Results.Trial.Situation(iT) = situation;
    Results.Trial.CUE_ARD(iT) = CUE_ARD;
    Results.Trial.CUE_MTL(iT) = CUE_MTL;
    Results.Trial.flipTime(iT) = flipTime;
    Results.Trial.FB_ARD(iT) = FB_ARD;
    Results.Trial.FB_MTL(iT) = FB_MTL;
    Results.Trial.fbTime(iT) = fbTime;
    Results.Trial.outcomeCode(iT) = outcome;
    Results.Session.Cycle = iT;
    if isStopTrial(iT)
        Results.Trial.STOP_ARD(iT) = STOP_ARD;
        Results.Trial.STOP_MTL(iT) = STOP_MTL;
        Results.Trial.stopTime(iT) = stopTime;
    else
        Results.Trial.STOP_ARD(iT) = 0;
        Results.Trial.STOP_MTL(iT) = 0;
        Results.Trial.stopTime(iT) = 0;
    end
    TrialLength = GetSecs - TrialLength;
    Results.Trial.TrialLength = TrialLength;
    Results.Trial.BackUp = backup;
    Results.Trial.BackUp_buff = backup_buff;
    save(results_fnm,'Results');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF TRIAL LOOP, SAVING, EXITING THE PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% End of experiment, saving the results
Screen('TextSize',win, 70);
DrawFormattedText(win, 'Vége a vizsgálatnak!', 'center', 'center', 0);          % Cue presentation
Screen(win,'Flip');
WaitSecs(2);
sca;    % Screen Close All
delser;
% keyboard;
% PDresplot('filename', results_fnm);
% keyboard