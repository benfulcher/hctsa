function s = gen(mode, len, frequency, samplerate, amplitude, startphase, yunit)

% signal/gen
%
% generate various kind of signals with time axis
%
% si = gen(mode, len, frequency, samplerate, amplitude, startphase, yunit)
%
% len is measured in seconds
% frequency, samplerate in Hertz
% amplitude is measured in the current yunit (default is dimensionless)
% startphase in radians
% yunit is an optional argument
%
% mode may be one of the following :
% sine              - sinusodial signal between [-amplitude +amplitude]
% unoise            - uniform distrubuted random numbers in interval [-amplitude +amplitude]
% gnoise            - gaussian random numbers with zero mean
% fnoise            - spectral flat noise with zero mean (DC value = 0)
% rect				- rectangular signal with values discret values 0 and amplitude
% symmrect          - rectangular signal with values discret values -amplitude and amplitude
% comb              - train of (one sample wide) spikes of height amplitude, otherwise signal is zero
% 
%
% cmerk Feb. 1998

if nargin < 1
	mode = 'sine';		% sine, uniform noise, gaussian noise, rect, puls ...
end
if nargin < 2
	len = 1;		% one second 
end
if nargin < 3
	frequency = 1000;	% 1000 Hz
end
if nargin < 4
	samplerate = 8000;	% Hertz
end
if nargin < 5
	amplitude = 1.0;	
end
if nargin < 6
	startphase = 0; 	% radian, 0..2pi
end
if nargin < 7
	yunit = unit;		% no dimension
end


len = ceil(len*samplerate);
xfirst = 0;
optinaltext = '';
switch mode
case	'sine'
	t = (2*pi*frequency/samplerate)*(0:len-1)' + startphase;
	tmp = amplitude * sin(t);
	optinaltext = [' with frequency : ' num2str(frequency) 'Hz'];
case	'unoise'
	tmp = (2* amplitude) * (rand(len,1)-0.5);
case	'gnoise'
	tmp = amplitude * gauss(len);
case	'fnoise'
	if mod(len, 2)
		N = (len-1)/2;
		r = exp(2*pi*i*rand(N,1));
		tmp = [0; r; conj(r(end:-1:1))];
	else							% even length 
		N = (len/2)-1;
		r = exp(2*pi*i*rand(N,1));
		tmp = [0; r; 1; conj(r(end:-1:1))];
	end
	tmp = (amplitude/2) * real(ifft(tmp));
case	'rect'
	tolerance = 1e-8;
	t = (2*pi*frequency/samplerate)*(0:len-1)' + startphase;
	tmp = sin(t);
	ind = find(tmp>tolerance);
	tmp(ind) = amplitude;
	ind = find(abs(tmp)<=tolerance);
	tmp(ind(1:2:end)) = amplitude;	% zeros values are set to amplitude zero amplitude zero ...
	tmp(ind(2:2:end)) = 0;
	ind = find(tmp<tolerance);
	tmp(ind) = 0;
	optinaltext = [' with frequency : ' num2str(frequency) 'Hz'];
case	'symmrect'
	tolerance = 1e-8;
	t = (2*pi*frequency/samplerate)*(0:len-1)' + startphase;
	tmp = sin(t);
	ind = find(tmp>tolerance);
	tmp(ind) = amplitude;
	ind = find(abs(tmp)<=tolerance);
	tmp(ind(1:2:end)) = amplitude;	% zeros values are set to amplitude zero amplitude zero ...
	tmp(ind(2:2:end)) = 0;
	ind = find(tmp<tolerance);
	tmp(ind) = -amplitude;
	optinaltext = [' with frequency : ' num2str(frequency) 'Hz'];
case	'comb'
	t = (pi*frequency/samplerate)*(0:len-1)' + startphase;	
	tmp = sin(t);
	zerocrossing = tmp(1:end-1) .* tmp(2:end);
	ind = find(zerocrossing<0);
	tmp = zeros(len,1);
	tmp(ind) = amplitude;
	if mod(startphase,2*pi) == 0
		tmp(1) = amplitude;
	end
	optinaltext = [' with frequency : ' num2str(frequency) 'Hz'];
case	'sinc'
	tolerance = 1e-10;
	t = (2* pi*frequency/samplerate)*(0:len-1)' + startphase;
	ind = find(abs(t)<tolerance);
	t(ind) = 1;
	tmp = amplitude * (sin(t)./t);
	tmp(ind) = amplitude;
	optinaltext = [' with frequency : ' num2str(frequency) 'Hz'];
otherwise	
	error(['mode ' mode ' not supported']);
end

s = signal(tmp, ['Generated ' mode ' signal' optinaltext]);
s = setaxis(s, 1,achse(unit('s'), xfirst, 1/samplerate));	
s = setyunit(s, yunit);
s = settype(s, 'scalar timeseries');
