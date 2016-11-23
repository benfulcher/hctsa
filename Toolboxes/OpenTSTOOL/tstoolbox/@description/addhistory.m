function d = addhistory(d, text)

%tstoolbox/@description/addhistory
%   adds text to current history list always the current time and date is
%   written into the first line
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,2);


if iscellstr(text)
	text{1} = [datestr(now) ' : ' text{1}];
else
	text = [datestr(now) ' : ' text];
end

d.history = TSTL_append(d.history, text);
