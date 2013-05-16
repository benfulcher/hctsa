function cout = surrogate1(cin)

%tstoolbox/@core/surrogate2
%   Syntax:
%     * cout = surrogate2(cin)
%
%   Input Arguments:
%     * cin - core object
%
%   create surrogate data for a scalar time series
%   see : James Theiler et al.'Using Surrogate Data to Detect Nonlinearity
%   in Time Series', APPENDIX : ALGORITHM II
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

if ndim(cin) > 1
	error('not a scalar time series')
end

x = data(cin);
N = dlens(cin,1);

[sx,index] = sort(x);		% sort time series
[dummy, rx] = sort(index);	% make ranked time series rx 

g = gauss(N);				% create gaussian random numbers
sg = sort(g);
y = sg(rx); 				% produce new time series out of origianl time series x,
							% but with gaussian amplitude distribution
							
y2 = data(surrogate1(core(y))); % use Theiler's Alg. I to make a surrogate of this gaussian time series
[dummy,index] = sort(y2);		
[dummy, ry2] = sort(index);	% make ranked time series ry2

x2 = sx(ry2);
cout = core(x2);




