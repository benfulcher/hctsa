function rs = realifft(data, N)

% inverse fft for fourier coefficents from real valued  original data
% needs the length of the original time series from which the fft was computed
% see also : realfft

data = data(:);
if mod(N,2) 	% odd length
	data = [data ; conj(data(end:-1:2))];  
else			% even length
	data = [data ; conj(data(end-1:-1:2))];  
end 

rs = real(ifft(data)); % remove small imaginary components introduced by rounding errors
