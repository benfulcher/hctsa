function filllistbox(datafiles,handles,fname)

global TSTOOLdatapath;

answer=get(handles.loadhandle,'UserData');
TSTOOLfilter=answer{1};
cwd=answer{2};


lhandle = handles.lboxhandle;

if ~isempty(datafiles)
  n=length(datafiles(:,1));
  m=length(datafiles(1,:));
else
  n=0;
  m=0;
end

String = get(lhandle, 'String');

String = {};	% clear existing contents in the listbox
activate=1;
for i=1:n
  for i1=1:m
    if ~isempty(char(datafiles(i,i1)))
%      disp(['nicht leer >' char(datafiles(i,i1)) '<']);
      last='';
      for i2=1:i1
	last=[last '      '];
      end	      
      
      if strcmp(char(datafiles(i,i1)),fname)
	activate=i;
      end
      last=[last char(datafiles(i,i1))];
    end
  end

  if isempty(String)
	String = {last};	% create a new cell array with one element 
  else
    String{end+1} = last;		% append name    
  end
end

%disp(['Name is >' fname '< at position ' num2str(activate)]);

set(lhandle, 'String', String); 	% make newly inserted name selected
set(lhandle, 'Value', activate);    
set(lhandle, 'UserData',datafiles);

%save(fullfile(TSTOOLdatapath,'tstool.mat'),'datafiles','TSTOOLfilter');
savesettings(handles);
