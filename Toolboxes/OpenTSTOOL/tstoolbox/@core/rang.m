function cout = rang(cin)

%tstoolbox/@core/rang
%   Syntax:
%     * cout = rang(cin)
%
%   Input Arguments:
%     * cin - core object
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

N = dlens(cin,1); 	% Laenge in der ersten Dimension (sollte auch die einzige sein (skalar))
ts = data(cin);		% Die Datenwerte des Input-Signals werden nach ts kopiert


[dummy, index] = sort(ts);
ts(index) = 1:N;			% reuse variable ts

cout = core(ts(:));
