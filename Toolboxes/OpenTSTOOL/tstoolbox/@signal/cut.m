function rs = cut(s, dim, start, stop)

%tstoolbox/@signal/cut
%   Syntax:
%     * rs = cut(s, dim, start, stop)
%
%   Input arguments:
%     * dim - dimension along which the signal is cutted
%     * start - position where to start the cut
%     * stop - position where to stop (optional)
%
%   Cut a part of the signal. If stop is ommited only the data at start is
%   cutted.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(3,4);

if (dim < 1) | (dim > ndim(s))
	error('Index dim exceeds matrix number of dimensions');
end

if (start < 1) | (start > dlens(s,dim) ) 
	error('Index start exceeds matrix dimensions');
end

if (nargin < 4)
	stop = start;
end

if (stop < 1) | (stop > dlens(s,dim)) | (stop < start) 
	error('Index stop exceeds matrix dimensions');
end


newlen = stop - start + 1;
points = data(s);
command = 'points = points('; 
for n=1:ndim(s)
	if n==dim
		command = [command 'start:stop,'];
	else
		command = [command ':,'];
	end
end
command = [command(1:end-1) ');'];
eval(command);

if newlen == 1
	if ndim(s) <= 2
		if dim==1
			points = points(:);
		end
	else
		points = squeeze(points);
	end
	rs = s;
	rs.core = core(points);
	rs.xaxes(dim) = [];
	%rs = removeaxis(rs, dim);
	
	a = getaxis(s, dim);
	a = cut(a, start,stop);
	lbl = label(a);
	v = spacing(a, 1); v = num2str(v(1));
	
	rs = addhistory(rs, ['(cut) Selected index ' num2str(start) ...
	' (' v ' ' lbl ') in dimension ' num2str(dim)]);	
else
	rs = s;
	rs.core = core(points);
	a = getaxis(rs, dim);
	a = cut(a, start,stop);
	lbl = label(a);
	v = spacing(a, stop-start+1);
	rs = setaxis(rs, dim, a);
	rs = addhistory(rs, ['(cut) Cut from index ' num2str(start) ' (' num2str(v(1)) ...
	' ' lbl ') to ' num2str(stop) ' (' num2str(v(end)) ' ' lbl ') in dimension ' num2str(dim)]);
end

rs = addcommandlines(rs, 's = cut(s', dim, start, stop);
