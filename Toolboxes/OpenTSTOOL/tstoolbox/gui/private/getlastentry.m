function name=getlastentry(datafiles,row)

if isempty(datafiles)
  name='';
else
  name=char(datafiles(row));
  for i=1:length(datafiles(1,:))
    if ~strcmp(datafiles(row,i),'')
      name=char(datafiles(row,i));
    end  
  end
end
  