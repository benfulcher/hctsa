function rs = amutual(s, maxtau, bins)

%tstoolbox/@signal/amutual
%   Syntax:
%     * amutual(s, maxtau, bins)
%
%   Input arguments:
%     * maxtau - maximal delay (should be much smaller than the lenght of
%       s) (optional)
%     * bins - number of bins used for histogram calculation (optional)
%
%   Auto mutual information function for real scalar signals, can be used
%   to determine a proper delay time for time-delay reconstruction. The
%   default value for maxtau is 25% of the input signal's length. The
%   default number of bins is 128.
%
%               [missing equation, see html/pdf documentation]
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,3);

if (ndim(s) > 1) | (~isreal(data(s)))
    error('Input signal must be a scalar time series');
end

N = dlens(s,1);

if nargin < 2
    maxtau = ceil(N/40);    % default : use 2.5 percent of input signals length      
end

if nargin < 3
    bins = 128;     
end

x = data(s);
[y,i] = sort(x);        % transformation to rang values
y(i) = (0:N-1)/N;		% and scaling to [0 1)

x = amutual(y, maxtau, bins);
rs = signal(core(x), s);	% special constructor calling syntax for working routines
a = getaxis(s, 1); 
a = setfirst(a, 0);
rs = setaxis(rs, 1, a);
rs = setyunit(rs, unit('Bit'));	
rs = addhistory(rs, ['Auto mutual information']);
rs = addcommandlines(rs, 's = amutual(s', maxtau, bins);
