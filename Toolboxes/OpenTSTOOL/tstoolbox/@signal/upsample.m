function rs = upsample(s, factor, method)

%tstoolbox/@signal/upsample
%   Syntax:
%     * rs = upsample(s, factor, method)
%
%   Input Arguments:
%     * method may be one of the following :
%          + 'fft'
%          + 'spline'
%          + 'akima'
%          + 'nearest'
%          + 'linear'
%          + 'cubic'
%     * s has be to sampled equidistantly for fft interpolation
%
%   Change sample rate of signal s by one-dimensional interpolation
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,3);

if nargin < 3
	method = 'fft'; 	% default upsampling method is fft
end

if ndim(s) > 1			
	error('Upsampling only for one-dimensional signals');
end

if factor < 1
	error('Upsampling factor must be greater one');
end

N = dlens(s,1); 		% old length of time-series
newN = ceil(factor*N);	% new length 
a = getaxis(s,1);

switch method
	case 'fft'
		switch  resolution(a)
			case 'linear'
				c = core(interpft(data(s), newN));
				rs = signal(c, s);				% special constructor syntax for class signal
				aa = setdelta(a, delta(a)*N/newN);
				rs = setaxis(rs,  1, aa);
				rs = addhistory(rs,  ['Upsampled using fft'] );
			otherwise
				error('Data values must sampled equidistantly for fft interpolation');
		end
	case {'spline', 'nearest', 'linear', 'cubic'}
				x = spacing(s, 1);
				y = data(s);
				xx = linspace(x(1), x(end), newN);		% compute new spacing values
				aa = achse(unit(a), xx(1), xx(2) - xx(1)); 	% new achse, cloned from old one 
				c = core(interp1(x(:), y(:), xx(:), method));
				rs = signal(c, s);
				rs = setaxis(rs, 1, aa);
				rs = addhistory(rs,  ['Upsampled using cubic splines'] );
	case 'akima'
				x = spacing(s, 1);
				y = data(s);
				xx = linspace(x(1), x(end), newN);		% compute new spacing values
				aa = achse(unit(a), xx(1), xx(2) - xx(1)); 	% new achse, cloned from old one 
				c = core(akimaspline(x(:), y(:), xx(:)));	% here the mex file 'akimaspline' is called 
				rs = signal(c, s);
				rs = setaxis(rs, 1, aa);
				rs = addhistory(rs,  ['Upsampled using akima splines'] );
	otherwise
		error(['Upsampling method ''' method ''' not supported']);	
end		

rs = addcommandlines(rs, 's = upsample(s', factor, method); 	% finally update commandlines 		

