function rs = aok(s, fensterlen, fftlen, vorschub, volume)

% signal/aok
%
% aok(s, fensterlen, fftlen, vorschub, volume)
%
% TFR (time-frequency representation) for a scalar signal s
% using an adaptive optimal kernel
%
% fensterlen - window length in samples
% fftlen -  fft length
% vorschub - moving step in samples
% volume - ?
% 
% view(aok(sine(2000, 1000, 8000))
%
% see also : spec2
% C.Merkwirth,U.Parlitz,W.Lauterborn  DPI Goettingen 1998
narginchk(1,5);
if nargin < 1, help(mfilename); return, end

if nargin < 2
      	fensterlen = 64;
        fftlen = 256;
        vorschub = 1;
		volume = 2;
elseif nargin < 3
        fftlen = 256;
        vorschub = 1;
		volume = 2;
elseif nargin < 4
        vorschub = 1;
		volume = 2;
elseif nargin < 5
		volume = 2;
end

c = aok(data(s), fensterlen, fftlen, vorschub, volume);
c = core(abs(c(:, ceil(fftlen/2):end)));
rs = signal(c, s);      % special constructor calling syntax for working routines
a = getaxis(s, 1);
rs = setaxis(rs, 2, achse(1/unit(a), 0, samplerate(a)/fftlen));
a = setdelta(a, vorschub*delta(a));
rs = setaxis(rs, 1, a);
rs = setplothint(rs, 'spectrogram');
rs = addhistory(rs, ['Calculated AOK TFR with window length ' num2str(fensterlen)]);
rs = addcommandlines(rs, 's = aok2(s', fensterlen, fftlen, vorschub, volume);


