function status = lclear(handles)

% tstool/lclear
% function status = lremove(handles)
% Loeschen aller Eintraege aus der Listbox und den UserData-Field der Figure, aber keine physik. Entfernung
% der Files

lhandle = handles.lboxhandle;

String = {};
set(lhandle, 'Value', 1);
set(lhandle, 'String', String);
setcurrentfile('', handles.lboxhandle, handles.currfilehandle); % jetzt ist kein File mehr der current working file
status = 0;



