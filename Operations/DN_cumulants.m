function out = DN_cumulants(y,n)
% Very simple function that uses the skewness and kurtosis functions in Matlab's Statistics Toolbox to calculate these higher order moments of input time series, y
% The appropriate higher order moment is specified with the second input, n
% Ben Fulcher, 2008

if nargin < 2 || isempty(n)
    error('DN_cumulants: you must specify a second input!')
end

switch n
	case 'skew1' % skewness
		out = skewness(y);
    case 'skew2' % corrects for bias
        out = skewness(y,0);
	case 'kurt1' % kurtosis
		out = kurtosis(y);
    case 'kurt2' % corrects for bias
        out = kurtosis(y,0);        
end


end