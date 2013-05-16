function rs = intspikeint(s, len)

%tstoolbox/@signal/intspikint
%   Syntax:
%     * rs = intspikeint(s)
%
%   Compute the interspike intervalls for a spiked scalar timeseries,
%   using transformation on ranked values.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

a = getaxis(s,1);
scalefac = delta(a);

N = dlens(cin,1); 	% Laenge in der ersten Dimension (sollte auch die einzige sein (skalar))
ts = data(cin);		% Die Datenwerte des Input-Signals werden nach ts kopiert

[dummy, index] = sort(ts);
ts(index) = 0:1:(N-1):1;			


c = core(ts(:));
rs = signal(c, s);				% special constructor calling syntax for working routines


rs = setyunit(rs, unit(a));
rs = setyname(rs, name(a));
rs = setaxis(rs, 1, achse);
rs = addhistory(rs,  ['Calculated interspike intervalls'] );
rs = addcommandlines(rs, 's = intspikeint(s', len);
