function d = settype(d, nt)

%tstoolbox/@description/settype
%   Syntax:
%     * d = settype(d, string)
%
%   Set a new type for signal
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

error(nargchk(2,2,nargin));

if isa(nt, 'char')
        d.type = nt;
else
        error('Wrong type of argument(s) given');
end


