function cout = surrogate1(cin)

%tstoolbox/@core/surrogate1
%   Syntax:
%     * cout = surrogate1(cin)
%
%   Input Arguments:
%     * cin - core object
%
%   create surrogate data for a scalar time series by randomizing phases
%   of fourier spectrum
%   see : James Theiler et al.'Using Surrogate Data to Detect Nonlinearity
%   in Time Series', APPENDIX : ALGORITHM I
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

if ndim(cin) > 1
	error('not a scalar time series')
end

y = data(cin);
N = dlens(cin,1);

z = realfft(y);

if mod(N,2)         % ungerade Länge
	nh = (N-1)/2;                  % N = 7  =>  nh = 3 
 	z(2:end) = z(2:end) .* exp(2*pi*i*rand(nh,1));
else
	nh = (N/2)-1;                  % N = 6  =>  nh = 2 
	z(2:end-1) = z(2:end-1) .* exp(2*pi*i*rand(nh,1));
end

cout = core(realifft(z,N)); 

