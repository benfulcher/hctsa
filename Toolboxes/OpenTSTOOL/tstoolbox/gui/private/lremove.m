function status = lremove(handles)

% tstool/lremove
% function status = lremove(handles)
% Loeschen des aktuellen Eintrags aus der Listbox und den UserData-Field der Figure, aber keine physik. Entfernung
% des Files

global TSTOOLdatapath;;



lhandle = handles.lboxhandle;
answer=get(handles.loadhandle,'UserData');
TSTOOLfilter=answer{1};
cwd=answer{2};
String = get(lhandle, 'String');
datafiles=get(lhandle,'UserData');

if ~isempty(String)
	Value =get(lhandle, 'Value');
	datafiles(Value,:)=[];
	String(Value) = []; 	% Dies ist der Syntax, um ein Element in einem Cell-Array zu loeschen
else
	Value = 1;
end

newsize = length(String);

if isempty(String)		% letzter Eintrag wurde geloescht
	set(lhandle, 'Value', 1);
	set(lhandle, 'String', String);
	status = 'no entry to remove';
elseif Value <= newsize
	set(lhandle, 'Value', Value);
	set(lhandle, 'String', String);
	String = get(lhandle, 'String');
	Value =get(lhandle, 'Value');	% aus Uebervorsicht wird  nochmal gecheckt
	filename = String{Value};		% so muss der Name wieder retrieved werden, nicht mittels ()
	status = 0;	
else
	Value = newsize;	% letzter Eintrag wird zum current working file
	set(lhandle, 'Value', Value);
	set(lhandle, 'String', String);
	String = get(lhandle, 'String');
	Value = get(lhandle, 'Value');	% aus Uebervorsicht wird  nochmal gecheckt
	filename = String{Value};		% so muss der Name wieder retrieved werden, nicht mittels ()
	status = 0;		
end

set(lhandle,'UserData',datafiles);

savesettings(handles);
%save(fullfile(TSTOOLdatapath,'tstool.mat'), 'datafiles','TSTOOLfilter');
