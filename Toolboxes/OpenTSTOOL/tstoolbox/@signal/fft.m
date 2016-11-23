function rs = fft(s)

%tstoolbox/@signal/fft
%   Syntax:
%     * f = fft(s)
%
%   Output arguments:
%     * f - n by 2 array, the first column contains the magnitudes, the
%       second one the phases.
%
%   Fourier transform of scalar signal s.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,1);

if ndim(s) > 1
	error('Only for scalar signals');
end

x = fftshift(fft(data(s)));
rs = signal(core([abs(x)/dlens(s,1) (abs(x) > 10 * sqrt(eps)) .* angle(x)]), s);
a = getaxis(rs, 1);
rs = setaxis(rs, 1, achse(unit(a)^(-1), -samplerate(a)/2, samplerate(a)/dlens(s,1)));
rs = setaxis(rs, 2, achse);
rs = setplothint(rs, 'subplotgraph');
rs = addhistory(rs, ['Fourier transform']);
rs = addcommandlines(rs, 's = fft(s');

