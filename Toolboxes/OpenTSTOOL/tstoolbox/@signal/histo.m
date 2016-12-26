function rs = histo(s, partitions)

%tstoolbox/@signal/histo
%   Syntax:
%     * histo(s, partitions)
%
%   Histogram function using equidistantly spaced partitions.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,2);

if (ndim(s) > 1) | (~isreal(data(s)))
	help(mfilename)
	return
end

epsilon = 1e-10;
N = dlens(s,1);

if (nargin < 2) 
	partitions = ceil(sqrt(N)/2);
end

points = data(s);
mi = min(points);
ma = max(points);

points = 1+floor((points-mi)/((ma-mi)/(partitions-epsilon)));
hist = sparse(points, ones(N,1), 1/N);

c = core(full(hist));
rs = signal(c, s);	% special constructor calling syntax for working routines
a = achse(yunit(s), mi, (ma - mi)/ (partitions-1));
rs = setaxis(rs, 1, a);
rs = setyunit(rs, unit);		% pdf values are scalars without unit
rs = addhistory(rs, 'Computed histogramm');
rs = setplothint(rs, 'bar');
rs = addcommandlines(rs, 's = histo(s', partitions);
