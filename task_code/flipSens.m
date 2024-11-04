function [ARD_time, MTL_time, flipTime] = flipSens(serial, win, flipInx)
%FLIPSENS Executing Flips for CUE, STOP, and FB presentation
%  Detecting photosensor timestamps
%  Returning ARD timestamps, MTL serial read time, MTL Fliptime
% FlipInx: 1 - Cue, 2 - Flipped, 3 - , 4 - Stop signal
flipTime = [];
global Sinx
global backup
global backup_buff
buttoncode = [];
terminator = 1;   % Error2 terminator if no more data coming from the serial port
termT = [];   % Error cycle terminator in case of no Flip detection for 0.2s
termT2 = [];
while isempty(buttoncode)
    [res, MTL_time, fliptime] = fgetl_nonblocking(serial, flipInx, win);
    if isempty(termT2)
        termT2 = GetSecs;
    end
    if ~isempty(fliptime)
        flipTime = fliptime;   % set it to output only when serial detection occured
    end
    
    if res   % One line of data from the serial port
        buttoncode = sscanf(res, '%d: %d');
        press_t = findstr(res, ' ');
        try
            ARD_time = str2num(res(press_t(2)+1:press_t(3)-1));
        catch
            buttoncode = [];
            press_t = [];
            disp('Flipsens chunk')
        end
    end
    
    if buttoncode ~= 4   % In case of not the sensor data are detected
        if isempty(flipTime)   % NOT FLIPPED
            if flipInx == 4 % In case it is a stop cycle
                disp('ERROR1: KEYPRESS BEFORE STOP FLIP')
                disp(res)
                Sinx.case = 1;
                Sinx.buttoncode(end+1) = buttoncode;
                Sinx.ARD_time(end+1) = ARD_time;
                Sinx.MTL_time(end+1) = MTL_time;
            end
        else
            if flipInx == 4
                disp('ERROR1.2: ALREADY FLIPPED')
                flipInx = 2;
            end
        end
        
        if flipInx == 2   % ALREADY FLIPPED
            disp('ERROR2: KEYPRESS AFTER STOP, BEFORE SENSOR DATA')
            disp(res)
            terminator = terminator + 1; % Counts backup fill
            Sinx.case = 2;
            Sinx.buttoncode(end+1) = buttoncode;
            Sinx.ARD_time(end+1) = ARD_time;
            Sinx.MTL_time(end+1) = MTL_time;
            termT = MTL_time;
            buttoncode = [];
            ARD_time = [];
            MTL_time = [];
            if ~isempty(find(backup_buff==13))
                currstr = backup_buff;
                buffer = currstr;
                endPos = find(buffer==13);
                buffsize = length(endPos);
                if buffsize > 1
                    if buffsize == 2
                        disp(['WARNING: TWO EVENT IN THE BACKUP BUFFER: ' num2str(buffsize)])
                    else
                        disp(['WARNING: MULTIPLE EVENT IN THE  BACKUP BUFFER: ' num2str(buffsize)])
                    end
                end
                for i = 1:buffsize
                    if backup(end-(buffsize-i),1) == 4
                        disp('FLIP TIMESTAMP FOUND')
                        buttoncode = backup(end-(buffsize-i), 1);
                        ARD_time = backup(end-(buffsize-i), 2);
                        MTL_time = backup(end-(buffsize-i), 6);
                    else
                        disp('BACKUP_BUFF LOAD TO SINX')
                        Sinx.buttoncode(end+1) = backup(end-(buffsize-i), 1);
                        Sinx.ARD_time(end+1) = backup(end-(buffsize-i), 2);
                        Sinx.MTL_time(end+1) = backup(end-(buffsize-i), 6);
                    end
                end
                
            else
                buttoncode = [];
                ARD_time = [];
                MTL_time = [];
                if terminator > 4 % More than 4 events in backup
                    disp('BUFFER FILLED WITH MULTIPLE KEYPRESSES, CONTINUE ANYWAY')
                    buttoncode = Sinx.buttoncode(1);
                    ARD_time = Sinx.ARD_time(1);
                    MTL_time = Sinx.MTL_time(1);
                    Sinx.flipTime = flipTime;
                    flipTime = NaN;
                end
            end
        else
            if flipInx~=4 % Delete presses during timeout
                disp('KEYPRESSES BEFORE CUE DETECTION')
                buttoncode = [];
                press_t = [];
                ARD_time = [];
                flipInx = 2;
            end
        end
    end
    if ~isempty(termT)
        if flipInx == 4
            if GetSecs - termT > 1.95
                disp('NO STOP FLIP DETECTED AFTER KEYPRESS, CONTINUE ANYWAY')
                buttoncode = Sinx.buttoncode(1);
                ARD_time = Sinx.ARD_time(1);
                MTL_time = Sinx.MTL_time(1);
                Sinx.flipTime = flipTime;
            end
        else
            if GetSecs - termT > 0.2
                disp('NO CUE FLIP DETECTED AFTER KEYPRESS, CONTINUE ANYWAY')
                buttoncode = Sinx.buttoncode(1);
                ARD_time = Sinx.ARD_time(1);
                MTL_time = Sinx.MTL_time(1);
                Sinx.flipTime = flipTime;
            end
        end
    end
    
    if flipInx == 4
        if GetSecs - termT2 > 1.95
            if fliptime == 0
                disp('Already forced to continue')
            else
                disp('NO STOP FLIP DETECTED, CONTINUE ANYWAY AFTER STOP')
                buttoncode = 4;
                ARD_time = 0;
                MTL_time = flipTime;
            end
        end
    else
        if GetSecs - termT2 > 0.2
            if fliptime == 0
                disp('Already forced to continue')
            else
                disp('NO CUE FLIP DETECTED, CONTINUE ANYWAY')
                buttoncode = 4;
                ARD_time = 0;
                MTL_time = flipTime;
            end
        end
    end
    flipInx = 2; % Already flipped
end
disp('DETECTED CUE: ')
disp(res)

