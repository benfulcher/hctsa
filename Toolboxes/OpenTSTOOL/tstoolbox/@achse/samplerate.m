function rate = samplerate(a)

%tstoolbox/@achse/samplerate
%   Syntax:
%     * rate = samplerate(a)
%
%   samplerate returns samplerate of achse object.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

if a.delta > 0
	rate = 1/a.delta;
else
	rate = 0;
end
