function rs = acp(s, tau, past, maxdelay, maxdim, nref)

%tstoolbox/@signal/acp
%   Syntax:
%     * acp(s, tau, past, maxdelay, maxdim, nref)
%
%   Input arguments:
%     * tau - proper delay time for s
%     * past - number of samples to exclude before and after reference
%       index (to avoid correlation effects)
%     * maxdelay - maximal delay (should be much smaller than the lenght
%       of s) (optional)
%     * maxdim - maximal dimension to use (optional)
%     * nref - number of reference points (optional)
%
%   Auto crossprediction function for real scalar signals for increasing
%   dimension. The default value for maxdelay is 25% of the input signal's
%   length. The default for maxdim is 8 and for nref it is 10% of the
%   input signal's length.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


narginchk(3,6);

if (ndim(s) > 1) | (~isreal(data(s)))
    error('Input signal must be a scalar time series');
end

N = dlens(s,1);

if nargin < 4
    maxdelay = max(16,ceil(N/40));    % default : use 2.5 percent of input signals length      
end
if nargin < 5
    maxdim = 8;     
end
if nargin < 6
    nref = max(100,ceil(N/10));   
end

NNR = max(2, ceil(N/1000));

x = zeros(maxdelay+1, maxdim);

ref = floor(rand(nref,1) * (N-maxdelay-maxdim*tau-1))+1;

for dim=1:maxdim
	E = data(embed(s, dim, tau));
	[a,b] = crossprediction(E(1:end-maxdelay-1,:), E, ref, maxdelay, NNR, past);
	x(:,dim) = b(:);
end

rs = signal(core(x), s);	% special constructor calling syntax for working routines
a = getaxis(s, 1); 
a = setfirst(a, 0);
rs = setaxis(rs, 1, a);
a2 = setname(achse(unit, 1, 1), 'Embedding dimension');
rs = setaxis(rs, 2, a2);
rs = setplothint(rs, 'multigraph');
rs = setyunit(rs, unit);	
rs = addhistory(rs, ['Auto cross prediction']);
rs = addcommandlines(rs, 's = acp(s', maxdelay, maxdim, nref);
