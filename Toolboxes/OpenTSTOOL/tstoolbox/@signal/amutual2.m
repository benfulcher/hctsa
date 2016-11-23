function rs = amutual2(s, len)

%tstoolbox/@signal/amutual2
%   Syntax:
%     * amutual2(s, len)
%
%   Input arguments:
%     * len - maximal lag
%
%   Auto mutual information (average) function for real scalar signals
%   using 128 equidistant partitions.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt
narginchk(2,2);

    if (ndim(s) > 1) | (~isreal(data(s)))
	help(mfilename)
	return
end

c = amutual2(s.core, len);
rs = signal(c, s);	% special constructor calling syntax for working routines
a = getaxis(s, 1); 
dl = delta(a);
a = setfirst(a, 0);
rs = setaxis(rs, 1, a);
rs = setyunit(rs, unit('Bit'));		% acf values are scalars without unit
rs = addhistory(rs, ['Auto mutual information of length ' num2str(len)]);
rs = addcommandlines(rs, 's = amutual(s', len);
