function cout = spec2(cin, fensterlen, fenster, vorschub)

%tstoolbox/@core/spec2
%   Syntax:
%     * cout = spec2(cin, fensterlen, fenster, vorschub)
%
%   Input Arguments:
%     * cin - core object
%     * fensterlen - window size
%     * fenster - type of window
%     * vorschub - moving step
%
%   spectrogramm of data using short time fft
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

win = window(fensterlen, fenster);
x = data(cin);

n = floor((length(x) - fensterlen)/vorschub);

tmp = zeros(n, fensterlen);

for i=1:n
	offset = (i-1)*vorschub;
	tmp(i, 1:fensterlen) = (x(offset+1:offset+fensterlen) .* win)'; 
end

tmp = abs(fft(tmp, [], 2));
cout = core(tmp(:,1:ceil(fensterlen/2)+1)/ceil(fensterlen/2));


