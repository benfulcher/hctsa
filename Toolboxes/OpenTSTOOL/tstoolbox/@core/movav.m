function cout = movav(cin, len)

%tstoolbox/@core/movav
%   Syntax:
%     * movav(cin,len)
%
%   Input Arguments:
%     * cin - core object
%     * len - average length
%
%   moving average
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

N = dlens(cin,1); 	% Laenge in der ersten Dimension (sollte auch die einzige sein (skalar))
ts = data(cin);		% Die Datenwerte des Input-Signals werden nach ts kopiert

av = zeros(N-len+1,1);

for i=1:N-len+1
	av(i) = mean(ts(i:i+len-1));
end

cout = core(av);
