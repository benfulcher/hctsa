function d = settype(d, nt)

%tstoolbox/@description/settype
%   Syntax:
%     * d = settype(d, string)
%
%   Set a new type for signal
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,2);

if isa(nt, 'char')
        d.type = nt;
else
        error('Wrong type of argument(s) given');
end


