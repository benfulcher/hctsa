function [status,datafiles] = siglevel(s,datafiles)

% print count of renamings.


str=history(s);
if ~isempty(datafiles)
  m=length(datafiles(:,1))+1;
else
  disp('There are no loaded files...');
  m=1;
end

n=0;
status=-1;
lasttoken='';
%disp(['add new datafiless ' name(s)]);
%disp(str);
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

if ~strcmp(lasttoken,name(s))
    if n      
      m=m+1;
      for i1=1:n
	datafiles(m,i1)=datafiles(m-1,i1);
      end
    end
    n=n+1;    
    datafiles(m,n)=cellstr(name(s));
end  
%disp(datafiles);

datafiles(cellfun('isempty',datafiles))=cellstr('');




return


