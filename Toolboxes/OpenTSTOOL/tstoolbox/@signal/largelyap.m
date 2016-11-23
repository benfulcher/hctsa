function [rs, s] = largelyap(s, n, stepsahead, past, nnr)

%tstoolbox/@signal/largelyap
%   Syntax:
%     * rs = largelyap(s, n, stepsahead, past, nnr)
%
%   Input arguments:
%     * n - number of randomly chosen reference points (-1 means: use all
%       points)
%     * stepsahead - maximal length of prediction in samples
%     * past - exclude
%     * nnr - number of nearest neighbours (optional)
%
%   Output arguments:
%     * rs -
%
%   Compute the largest lyapunov exponent of a time-delay reconstructed
%   timeseries s, using formula (1.5. of Nonlinear Time-Series Analysis,
%   Ulrich Parlitz 1998 ).
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(4,5);

if nargin < 4
	help(mfilename)
	return;
end

if nargin < 5
	nnr = 1;
end

[N,dim] = size(data(s));

ref = randref(1, N-stepsahead, n);

try 
	atria = optparams(s, 1);
	names = fieldnames(atria);
catch
	atria = nn_prepare(data(s), 'euclidian');
	s = setoptparams(s, 1, atria);
end

a = getaxis(s,1);
rs = signal(core(largelyap(atria, data(s), ref, stepsahead, nnr, past)), s);

a = setfirst(a,0);	
rs = setaxis(rs, 1, a);
rs = addhistory(rs,  ['Computed largest lyapunov exponent']);
rs = addcommandlines(rs, 's = largelyap(s', n, stepsahead, past, nnr);
rs = setyname(rs, 'p');
rs = setlabel(rs, 'Prediction error');

end