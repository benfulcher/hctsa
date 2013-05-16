function display(d)

%tstoolbox/@description/display
%   description/display
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt



disp(['  Name : ' d.name]);
disp(['  Type : ' d.type]);
disp(' ')
disp(['  Attributes of data values : '])
disp(['  ' label(d.yunit) ' | '  d.yname])
disp(' ')
disp('  Comment : ')
display(d.comment);
disp(' ')
disp('  History : ')
display(d.history);





