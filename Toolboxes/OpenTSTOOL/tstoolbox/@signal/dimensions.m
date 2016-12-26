function [rs, rs2, rs3] = dimensions(s, bins)

%tstoolbox/@signal/dimensions
%   Syntax:
%     * [bc,in,co] = dimensions(s, bins)
%
%   Input arguments:
%     * s - data points (row vectors)
%     * bins - maximal number of partition per axis, default is 100
%
%   Output arguments:
%     * bc - scaling of boxes with partititon sizes (log[2]-log[2])
%     * in - scaling of information with partititon sizes (log[2]-log[2])
%     * co - scaling of correlation with partititon sizes (log[2]-log[2])
%
%   Compute boxcounting, information and correlation dimension of a
%   time-delay reconstructed timeseries s for dimensions from 1 to D,
%   where D is the dimension of the input vectors using boxcounting
%   approach.
%
%   Scale data to be within 0 and 1. Give a sortiment of (integer)
%   partitionsizes with almost exponential behaviour.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,2);

if ndim(s) ~= 2
    error('Signal must contain vector data');    
end

if nargin<2
    bins = 100;    
end

points = data(s);
[N,dim] = size(points);

% scale data to be within 0 and 1
points = points - min(min(points));
points = points / max(max(points));

% give a sortiment of (integer) partitionsizes with almost exponential behaviour 
par = [2 3 4 5 6 7 8 10 12 14 16 20 23 27 32 39 46 54 64 77 91 108 128 153 182 216 256 ...
      305 363 431 512 609 725 862 1024 1218 1449 1723 2048 2436 2897 3445 4096 4871 ...
	  5793 6889 8192 9742 11586 13778 16384 19484 23171 27555 32768 38968 46341 55109 65536];

partitions = par(find(par<=bins));	% use no sizes greater than bins

[c,d,e] = boxcount(points, partitions);

c = [zeros(1,dim) ; c];             % add zeros for partition size 1
d = [zeros(1,dim) ; d];             % add zeros for partition size 1
e = [zeros(1,dim) ; e];             % add zeros for partition size 1
partitions = [1 ; partitions(:)];

a1 = achse(-log2(partitions));     		% create axis with arbitrary spacing
a1 = setname(a1, 'ld r');

a2 = setname(achse(unit, 1, 1), 'Embedding dimension');

rs = signal(core(c), s);	
rs = setaxis(rs, 1, a1);
rs = setaxis(rs, 2, a2);
rs = setplothint(rs, 'multigraph');
rs = addhistory(rs,  ['Computed boxcounting dimension']);
rs = addcommandlines(rs, 's = boxdim(s', bins);
rs = setyname(rs, 'ld N(r)');
rs = setlabel(rs, 'Scaling of D0');

rs2 = signal(core(d), s);	
rs2 = setaxis(rs2, 1, a1);
rs2 = setaxis(rs2, 2, a2);
rs2 = setplothint(rs2, 'multigraph');
rs2 = addhistory(rs2,  ['Computed information dimension']);
rs2 = addcommandlines(rs2, 's = infodim(s', bins);
rs2 = setyname(rs2, 'I(r)');
rs2 = setlabel(rs2, 'Scaling of D1');

rs3 = signal(core(e), s);	
rs3 = setaxis(rs3, 1, a1);
rs3 = setaxis(rs3, 2, a2);
rs3 = setplothint(rs3, 'multigraph');
rs3 = addhistory(rs3,  ['Computed correlation dimension']);
rs3 = addcommandlines(rs3, 's = corrdim(s', bins);
rs3 = setyname(rs3, 'ld C(r)');
rs3 = setlabel(rs3, 'Scaling of D2');

