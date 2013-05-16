function cout = acf(cin, m)

%tstoolbox/@core/acf
%   Syntax:
%     * acf(cin, m)
%
%   Input Arguments:
%     * cin core object
%     * m fft-length
%
%   acf calculates the autocorrelation function of cin via fft of length
%   m.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


x = data(cin);


if mod(m,2)
	m = m-1;%<-- make m even, because m = FFT length
	x(end) = [];
end
m2 = m/2;
Lx = length(x);
mu = (-1).^(0:m-1)';
ak = zeros(m,1);
z  = zeros(m,1);
n1 = 1;
while( n1 <= Lx )
   n2 = min( n1+m2-1, Lx );
   zi = fft(x(n1:n2), m);
   ak = ak + zi.*conj( zi + mu.*z );
   z  = zi;
   n1 = n1 + m2;
end;
ak = ifft(ak);
if( ~any(imag(x)) ),
   ak = real(ak);  end
ak = ak(1:m2+1);
cout = core(ak/ak(1));
	
