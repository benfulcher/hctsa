function [status,datafiles] = writelevel(filename,datafiles)

% print count of renamings.
%datafiles=cell(10,10);

s=signal(filename);
str=history(s);
if ~isempty(datafiles)
  m=length(datafiles(:,1))+1;
else
  m=1;
end
n=0;
status=-1;
lasttoken='';
for i=1:length(str(:,1))
  [token,rem]=strtok(str(i,:),'>');
  rem=rem(2:length(rem));
  [token,rem]=strtok(rem,'<');
  if token & ~strcmp(lasttoken,token)
    lasttoken=token;
    if n      
      m=m+1;
      for i1=1:n
	datafiles(m,i1)=datafiles(m-1,i1);
      end
    end
    n=n+1;    
    datafiles(m,n)=cellstr(token);
    status=0;
  end
end


[path,fname,ext,ver] = fileparts(filename);

if ~strcmp(lasttoken,fname)
    if n      
      m=m+1;
      for i1=1:n
	datafiles(m,i1)=datafiles(m-1,i1);
      end
    end
    n=n+1;    
    datafiles(m,n)=cellstr(fname);
end  

datafiles(cellfun('isempty',datafiles))=cellstr('');




return;


