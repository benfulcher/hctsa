function rs = spec2(s, fensterlen, fenster, vorschub)

%tstoolbox/@signal/spec2
%   Syntax:
%     * rs = spec2(s)
%
%   Input Arguments:
%     * fensterlen - size of window (optional)
%     * fenster - window type (optional)
%     * vorschub - shift in samples (optional)
%
%   spectrogramm of signal s using short time fft
%
%   Examples:
%view(spec2(sine(10000, 1000, 8000), 512, 'Hanning'))
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,4);

if nargin < 2
	if dlens(s,1) > 200
		fensterlen = 128;
	else
		fensterlen = ceil(dlens(s,1)/4);
	end
	fenster = 'Hamming';
	vorschub = fensterlen/4;
elseif nargin < 3
	fenster = 'Hamming';
	vorschub = fensterlen/4;
elseif nargin < 4
	vorschub = fensterlen/4;
end

win = window(fensterlen, fenster);
x = data(s);
size_win = size (win);
size_x = size (x);
if size_win(2) ~= size_x(2)
  errordlg (['for windowtype ' fenster ' signal must have ' ...
	     num2str(size_win(2)) ' columns'], 'wrong dimensions')
  rs = -1;
  return
end

c = spec2(s.core, fensterlen, fenster, vorschub);
rs = signal(c, s);	% special constructor calling syntax for working routines
a = getaxis(s, 1);
rs = setaxis(rs, 2, achse(1/unit(a), 0, samplerate(a)/fensterlen));
a = setdelta(a, vorschub*delta(a));
rs = setaxis(rs, 1, a);
rs = setplothint(rs, 'spectrogram');
rs = addhistory(rs, ['Calculated spectrogram with window length ' num2str(fensterlen)]);
rs = addcommandlines(rs, 's = spec2(s', fensterlen, fenster, vorschub);



