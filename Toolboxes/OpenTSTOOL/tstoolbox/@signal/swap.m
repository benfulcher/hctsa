function rs = swap(s, dim1, dim2)

%tstoolbox/@signal/swap
%   Syntax:
%     * rs = swap(s) (exchange dimension 1 and dimension 2)
%     * rs = swap(s, dim1, dim2)
%
%   Exchange signal's dimensions (and axes)
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,3);

if nargin == 2
	error('Need one or three arguments');
end

if nargin < 3
	dim1 = 1;
	dim2 = 2;
end

if (dim1 < 1) | (dim1 > ndim(s))
	error(['dim1 must be between 1 and' num2str(ndim(s))]);
end

if (dim2 < 1) | (dim2 > ndim(s))
	error(['dim1 must be between 1 and' num2str(ndim(s))]);
end

if dim1 == dim2
	rs = s;
	return;
end

order = 1:ndim(s);

order(dim1) = dim2;
order(dim2) = dim1;

rs = s;
rs.core = core(permute(data(s), order));
a1 = getaxis(s, dim1);
a2 = getaxis(s, dim2);
rs = setaxis(rs, dim1, a2);
rs = setaxis(rs, dim2, a1);

rs = addhistory(rs, ['Swapped dimension ' num2str(dim1) ' with dimension ' num2str(dim2)]);
rs = addcommandlines(rs, 's = swap(s', dim1, dim2);
