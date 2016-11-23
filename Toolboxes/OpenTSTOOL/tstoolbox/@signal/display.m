function display(s)

%tstoolbox/@signal/display
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt
narginchk(1,1);

disp(['  ' inputname(1) ' = signal object '])
disp(' ');
display(s.core)
for i=1:ndim(s)
	a = getaxis(s, i);
	disp(['  X-Axis ' num2str(i) ' : ' label(a) ' | ' name(a)]);
end
disp(' ');
display(s.description)

