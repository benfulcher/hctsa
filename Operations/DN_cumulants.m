function out = DN_cumulants(y,whatcum)
% Very simple function that uses the skewness and kurtosis functions in 
% Matlab's Statistics Toolbox to calculate these higher order moments of input time series, y
% The appropriate higher order moment is specified with the second input, n
% Ben Fulcher, 2008

if nargin < 2 || isempty(whatcum)
    whatcum = 'skew1'; % do skewness by default
end

switch whatcum
case 'skew1' % skewness
	out = skewness(y);
case 'skew2' % corrects for bias
    out = skewness(y,0);
case 'kurt1' % kurtosis
	out = kurtosis(y);
case 'kurt2' % corrects for bias
    out = kurtosis(y,0);        
otherwise
    error('Unknown cumulant specified: should be ''skew1'', ''skew2'', ''kurt1'', or ''kur2''')
end

end