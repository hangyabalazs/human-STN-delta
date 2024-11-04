function [buttoncode, ARD_time, MTL_time] = ser_read(serial, duration, Start, win, flipInx)
%SER_READ Keypress detection with FGETL_NONBLOCKING
%  Detecting buttonpress timestamps
%  Returning pressed button code,  ARD timestamps, MTL serial read time
global cue
global RectSize
global backup_buff
buttoncode = [];
res = [];
if duration == 2
    sTime = Inf;
else
    sTime = duration-(duration-0.1);
end
if Start==0 || isempty(Start)
    disp('No starttime inhereted')
    Start = GetSecs; % Fake timestamp forcing script to continue
end
while isempty(buttoncode) && GetSecs-Start < duration
    [res, MTL_time, ~] = fgetl_nonblocking(serial, flipInx, win);
    if res
        buttoncode = sscanf(res, '%d: %d');
        press_t = findstr(res, ' ');
        ARD_time = str2num(res(press_t(2)+1:press_t(3)-1));
    end
    
    if buttoncode == 4
        buttoncode = [];
        press_t = [];
        ARD_time = [];
        disp('FLIP DETECTED IN SER_READ. DELETE IT')
    end
    
    if GetSecs-Start > sTime
        sTime = Inf;
        Screen('TextSize',win, 300);
        DrawFormattedText(win, cue, 'center', 'center', 0);          % Cue presentation
        Screen('Fillrect', win, [255 255 255], [0 0 RectSize]);   % white box before stop signal for the sensor
        Screen('Flip', win);   % present to the screen
        disp('WHITE BOX')
    end
    
    if isempty(buttoncode)
        MTL_time = [];
        if ~isempty(find(backup_buff==13))
            disp('LOAD BACKUP_BUFF IN SER_READ')
            backup_buff = char(backup_buff);
            disp(backup_buff)
            buttoncode = sscanf(backup_buff, '%d: %d');
            press_t = findstr(backup_buff, ' ');
            while buttoncode == 4
                backup_buff = backup_buff(press_t(4)+1:end);
                disp(backup_buff)
                buttoncode = sscanf(backup_buff, '%d: %d');
                press_t = findstr(backup_buff, ' ');
            end
            if isempty(press_t)
                buttoncode = [];
                backup_buff = [];
                disp('Some piece in backup_buff. Erase')
            else
                ARD_time = str2num(backup_buff(press_t(2)+1:press_t(3)-1));
                press_t = backup_buff(press_t(4)+1:end);
                backup_buff = [];
            end
        end
    end
    
    
    
end
buttoncode = buttoncode + 1;   % set inputnumbers to keycodes

if isempty(buttoncode)   % In case of keypress withold
    ARD_time = [];
    MTL_time = [];
end

disp('DETECTED KEY: ')
if res
    disp(res)
elseif buttoncode
    disp([num2str(buttoncode-1) ' at ' num2str(ARD_time)  ' press ' num2str(press_t)])
elseif isempty(buttoncode)
    disp('EMPTY')
end