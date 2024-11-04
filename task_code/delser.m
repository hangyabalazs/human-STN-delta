function delser

delete(instrfindall())
ser = serial('COM3', 'Baudrate', 115200);
fopen(ser);
fclose(ser);   % close and delete serial comport
delete(ser);
clear ser;