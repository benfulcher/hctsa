function status = lsort(handles)

% tstool/lsort
% sorts the entries in the listbox with handle handles.lboxhandle
% and removes duplicate entries

if nargin < 1, help(mfilename), return, end 

lboxhandle = handles.lboxhandle;

Value = get(lboxhandle, 'Value');
String = get(lboxhandle, 'String');
if ~isempty(String) 
	[String, index] = cellsort(String);
	String = eliminatedupl(String);
	set(lboxhandle, 'String', String);
	set(lboxhandle, 'Value', 1);	
	setcurrentfile(String{1}, handles.lboxhandle, handles.currfilehandle);
	status = 0; 	% OK
else
	status = 'no files to sort';
end		
