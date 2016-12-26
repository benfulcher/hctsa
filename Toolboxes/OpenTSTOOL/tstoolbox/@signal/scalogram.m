function rs = scalogram(s, scalemin, scalemax, scalestep, mlen)

%tstoolbox/@signal/scalogram
%   Syntax:
%     * rs = scalogram(s) => scalemin=0.1
%     * rs = scalogram(s, scalemin) => scalemax=1
%     * rs = scalogram(s, scalemin, scalemax) => scalestep=0.1
%     * rs = scalogram(s, scalemin, scalemax, scalestep) => mlen=10
%     * rs = scalogram(s, scalemin, scalemax, scalestep, mlen)
%
%   Scalogram of signal s using morlet wavelet. See also: spec2.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,5);

if nargin<2
	scalemin = 0.1;
end
if nargin<3
	scalemax = 1;
end
if nargin<4
	scalestep = 0.1;
end
if nargin<6
	mlen = 10;
end

c = scalogram(s.core, scalemin, scalemax, scalestep, mlen);
rs = signal(c, s);	% special constructor calling syntax for working routines
a = achse(unit, scalemin, scalestep);
rs = setaxis(rs, 2, a);
rs = setplothint(rs, 'spectrogram');
rs = addhistory(rs, ['Calculated scalogram']);
rs = addcommandlines(rs, 's = scalogram(s', scalemin, scalemax, scalestep);
