function cout = rms(cin)

%tstoolbox/@core/rms
%   Syntax:
%     * cout = rms(cin)
%
%   Input Arguments:
%     * cin - core object
%
%   compute root mean square value of each column of c1
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

L = dlens(cin,1);
tmp = data(cin);
tmp = sqrt(sum(tmp.*conj(tmp),1)/L);

if (dlens(cin) > 1)
	cout = core(shiftdim(tmp,1));	% remove first dimension
else
	cout = core(tmp);
end
