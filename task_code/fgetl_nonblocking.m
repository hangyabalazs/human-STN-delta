function [ retval, MTL_time, flipTime ] = fgetl_nonblocking(serialport, flipInx, win)
%FGETL_NONBLOCKING Nonblocking read of serial data, returning entire lines
%  Collect characters from serial port, return empty array or a line of
%  string if line terminator (CR or LF) was hit
global backup
global backup_buff
persistent buff

if flipInx == 4 % Stop presentation
    if serialport.BytesAvailable > 0   % checks if bytes are already available in the buffer before Flip
        disp('STOP PREFLIP BUFFER CORRECTION')
        flipInx = 2;
    elseif ~isempty(find(buff==13))
        disp('STOP PREFLIP BUFFER ALREADY FILLED')
        flipInx = 5;
    end
end

switch flipInx
    case 1   % CUE
        if serialport.BytesAvailable > 0   % checks if bytes are already available in the buffer before Flip
            currstr = fread(serialport, serialport.BytesAvailable, 'char')';
            buff = [buff, currstr];
        end
        buff = [];
        backup_buff = [];
        Screen('Flip', win);   % CUE PRESENTATION TO THE SCREEN
        flipTime = GetSecs;   % Timestamp of the cue presentation
        disp('CUE FLIP')
    case 2   % SER READ
        flipTime = [];   % Only serial read, no cue presentation
    case 3   % FEEDBACK
        if serialport.BytesAvailable > 0   % checks if bytes are already available in the buffer before Flip
            currstr = fread(serialport, serialport.BytesAvailable, 'char')';
            buff = [buff, currstr];
        end
        buff = [];
        backup_buff = [];
        Screen('Flip', win);   % FEEDBACK PRESENTATION TO THE SCREEN
        flipTime = GetSecs;   % Timestamp of the feedback presentation
        disp('FEEDBACK')
    case 4   % STOP
        Screen('Flip', win);   % STOP PRESENTATION TO THE SCREEN
        flipTime = GetSecs;   % Timestamp of the stop presentation
        disp('STOP FLIP')
    case 5
        flipTime = [];
        if ~isempty(find(backup_buff==13))
            disp('LOSING DATA (5)')
            disp(char(buff))
            disp(char(backup_buff))
            if ~strcmp(char(buff), char(backup_buff))
                if ~isempty(find(buff==13))
                    backup_buff = buff;
                    disp('LOAD TO BACKUP_BUFF')
                    disp(char(backup_buff))
                else
                    disp('BACKUP_BUFF CHUNK DETECTED')
                    disp(backup_buff)
                end
            else
                disp('BACKUP_BUFF AND BUFF ARE THE SAME')
                if serialport.BytesAvailable > 0
                    backup_buff = [];
                end
            end
        end
end

retval = [];

if flipInx == 5
    MTL_time = GetSecs;
else
    MTL_time = [];
end

if serialport.BytesAvailable > 0
    MTL_time = GetSecs;
    disp('BUFFER CHECK: ')
    % compose a single line from the chunks read from serial port
    currstr = fread(serialport, serialport.BytesAvailable, 'char')';
    buffer = currstr;
    endPos = find(buffer==13);
    buffsize = length(endPos);
    if buffsize > 1
        if buffsize == 2
            disp(['WARNING: TWO EVENT IN THE BUFFER: ' num2str(buffsize)])
        else
            disp(['WARNING: MULTIPLE EVENT IN THE BUFFER: ' num2str(buffsize)])
        end
        buffer = char(buffer);
        disp(char(currstr))
        if strcmp(buffer(1:endPos(1)), buffer(endPos(1)+2:endPos(2)))
            first_CRLF = find(currstr==13 | currstr==10, 1);
            if first_CRLF
                currstr = buff(1:first_CRLF);
                disp('DELETE DUPLICATION. REMAINING: ')
                disp(char(currstr))
            end
        end
    end
    
    buff = [buff, currstr];
    
    buffer = buff;
    endPos = find(buffer==13);
    buffsize = length(endPos);
    
    
    for i = 1:buffsize
        if i < 2
        endLine = endPos(1)-1;
        else
            endLine = endPos(i)-(endPos(i-1)+2);
        end
        bLine = char(buffer(endPos(i)-endLine:endPos(i)));
        bData = sscanf(bLine, '%d %*s %d %*s %d');
        if length(bData) < 3
            disp('Data chunk in buffer')
            disp(bData)
        else
            iB = size(backup, 1)+1;
            backup(iB, 1) = bData(1);
            backup(iB, 2) = bData(2);
            backup(iB, 3) = bData(3);
            backup(iB, 4) = 1;
            backup(iB, 5) = flipInx;
            backup(iB, 6) = MTL_time;
            disp([num2str(bData(1)) ' at ' num2str(bData(2))  ' press ' num2str(bData(3))])
        end
    end
    
    if length(buff) > 100
        fprintf('Dropping buffer, length %d already over limits :', length(buff))
        disp(char(buff))
        buff = [];
    end
    
    % if there is carriage return or new line (CR == ASCII 13 or LF == ASCII 10) then
    % process until that
    
    first_CRLF = find(buff==13 | buff==10, 1);
    
    if first_CRLF
        buff;
        % return with string before line terminator,
        % keep buffer part following the terminator character
        retval = char(buff(1:first_CRLF));
        if buff(first_CRLF) == 13
            buff = buff(first_CRLF+1:end);
            if buff & buff(1) == 10
                % line end was a CR+LF (13+10), skip both
                buff = buff(2:end);
            end
        else
            % line end was either a CR or an LF
            buff = buff(first_CRLF+1:end);
        end
    end
    if ~isempty(find(buff==13))
        
        disp('REMAINING PIECE IN BUFFER')
        disp(char(buff))
        disp('RETVAL')
        disp(retval)
        if ~strcmp(retval, char(buff))
            backup_buff = buff;
            disp('LOAD TO BACKUP_BUFF2')
            disp(char(backup_buff))
        end
    end
end
end
