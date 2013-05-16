function [currentfile]=getcurrentfile(dir,name)

global TSTOOLfilter;
%  return the filename for the selected signal

if ~isempty(name)
  while (strcmp(name(1),' '))
    name=name(2:length(name));
  end
  
  currentfile=[fullfile(dir,name) '.sgo'];
else
  currentfile='';
end


return;

 
      
    

