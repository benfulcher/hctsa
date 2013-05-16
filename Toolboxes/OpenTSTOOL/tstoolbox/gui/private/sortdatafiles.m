function datafiles=sortdatafiles(datafiles)

% das ist ein Test
%

newdatafiles={};
new_n=0;
n=length(datafiles(:,1));
m=length(datafiles(1,:));

for i=1:n
  equal_line=0;
  for i1=1:i-1
    equal=1;
    for i2=1:m
      if ~strcmp(char(datafiles(i,i2)),char(datafiles(i1,i2)))
	equal=0;
      end	
    end
%    if equal	
%      disp(['Zeile ' num2str(i) ' und Zeile ' num2str(i1) ' sind' ...
%	      ' gleich!']);
%    end
    if equal
      equal_line=1;
    end
  end
  if ~equal_line
    new_n=new_n+1;
    for i2=1:m
      newdatafiles(new_n,i2)=datafiles(i,i2);
    end    
  end 
end

datafiles=newdatafiles;
newdatafiles={''};
neworder=[];
new_n=1;
for i=1:m                            % Anzahl der Elemente  
  lines=testelements(datafiles,new_n,i);
%  disp(' ');
%  disp([num2str(i) ' elements in lines ' num2str(lines)]);
  if i==1
    neworder=lines;
  else
    i2=0;
    while i2<length(neworder)      	
      i2=i2+1;
      suborder=[];
%      disp(' ');
      for i1=1:length(lines)
%	disp(['try to match line ' num2str(lines(i1)) ' and line ' num2str(neworder(i2))]);
	is_equal=1;
	for i3=1:i-1
	  if ~strcmp(char(datafiles(lines(i1),i3)),char(datafiles(neworder(i2),i3)))
	    is_equal=0;
%	    disp(['Theres a mismatch in element ' num2str(i3)]);
	  end
	end	
	if is_equal
%	  disp(['line ' num2str(neworder(i2)) ' correspond with line ' num2str(lines(i1))]);
	  suborder=[suborder lines(i1)];
	end	      
      end
      if i2<length(neworder)
	neworder=[neworder(1:i2) suborder neworder(i2+1: ...
						 length(neworder))];
      else
	neworder=[neworder(1:i2) suborder];
      end      
%      disp(['i2=' num2str(i2) ' New order ' num2str(neworder)]);
      i2=i2+length(suborder);	
    end
  end  
end

newdatafiles={};

for i=1:length(neworder)
  for i1=1:m
    newdatafiles(i,i1)=datafiles(neworder(i),i1);
  end
end  

datafiles=newdatafiles;
%disp(datafiles);

return

function [lines]=testelements(datafiles,n,m)
lines=[];
%disp(['finding row with ' num2str(m) ' elements:']);
for i=n:length(datafiles(:,1))
%  disp(['testing row ' num2str(i) '.']);
  number=i;
  for i1=1:m
    if isempty(char(datafiles(i,i1)))
      number=0;
%      disp(['element ' num2str(i1) ' is empty!']);
    end   
  end
  for i1=m+1:length(datafiles(1,:))
    if ~isempty(char(datafiles(i,i1)))
      number=0;
%      disp(['element ' num2str(i1) ' is not empty!']);
    end
  end  
  if number 
    lines=[lines number];
  end
end
if lines
%  disp('found it!');
end
return;



