function rs = realfft(data)

% fft for real valued data, returns only the necessary 
% fourier coefficients

N = length(data);

rs = fft(data(:));

if mod(N,2) 	% odd length
	nh = (N+1)/2;  			% N = 5  =>  nh = 3 
	rs = rs(1:nh);				% coefficient for hightest freq. is complex
else			% even length
	nh = (N/2)+1;  			% N = 6  =>  nh = 4 
	rs = rs(1:nh);				% coefficient for hightest freq. is real
end 
