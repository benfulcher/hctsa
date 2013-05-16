function cout = scalogram(cin, smin, smax, sstep, tim)

%tstoolbox/@core/scalogram
%   Syntax:
%     * cout = scalogram(cin, smin, smax, sstep, tim)
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

x = data(cin);
lx = dlens(cin,1);

s = smin:sstep:smax;
sc = zeros(lx, length(s));

index = 1;
for scp=smin:sstep:smax
	[w,tw]=morletw(8 * tim,scp,0,tim); 
	lw=length(w);

	 % The Morlet wavelet should be sampled at least at a 2 Hz rate with
	 % seconds as time unit, or 2 kHz rate while using milliseconds... If
	 % the scale is too low, the wavelet may result too much compressed
	 % and wavelet properties are lost due to undersampling.
	 % If this happens, the user should increase the lowest scale or 
	 % change the time units to a greater value.

	if lw<4*tim+1  % slowest allowed sampling rate
		warning(['scale ' num2str(scp) ' may be too low. Increase smin']);
	end

	 % If the scale is too high, the wavelet may result excesively dilated,
	 % and filtering operations may be too slow. This does not affect to
	 % the wavelet properties, but to the users patience...
	 % If this happens, decrease the highest scale or choose smaller 
	 % time units.
	 % Note that the higher is the scale, the larger is the wavelet.

	if lw>25000
		fprintf('\nWARNING: Scale (%f) may be too high. May take a long time ...\n',scp);
	end

	del=floor((lw-1)/2);
	w1=conv(w,x);
	y1=w1(del+1:del+lx);
	sc(:,index) = y1(:);
	index = index + 1;
end

cout = core(abs(sc));


function [w,t]=morletw(fs,a,b,ki)

%MORLETW  Calculate the Morlet Wavelet.
%
%         W = MORLETW (FS,SC) calculates the Morlet wavelet
%         at scale SC sampled at rate FS.
%
%         W = MORLETW (FS,SC,SHF) calculates the Morlet wavelet
%         at scale SC, sampled at rate FS and delayed SHF seconds.
%
%         W = MORLETW (FS,SC,SHF,TIM) calculates the Morlet wavelet
%         at scale SC, sampled at rate FS and delayed SHF seconds,
%         changing the time interval to -TIM:TIM seconds.
%


if ki<5
	warning('Small time interval, wavelet ends may be too abrupt');
end
	
% We find the number of points according to the sampling rate 
% and the time space.

n=ki*2*fs+1;
res=(2*ki)/n;             
N=(2*ki)*a/res;            % Number of points for the scaled wavelet

t=linspace(-ki*a,ki*a,N);  % Time domain for the scaled wavelet

wo=2*pi;                 % wo=2*pi makes the 1st scale have

w=sqrt(1/a)*exp(j*wo*(t-b)/a).*exp(-((t-b)/a).^2/2);
