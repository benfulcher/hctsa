% DN_Cumulants
% 
% Very simple function that uses the skewness and kurtosis functions in 
% Matlab's Statistics Toolbox to calculate these higher order moments of input time series, y
% 
% INPUTS:
% 
% y, the input time series
% 
% whatcum, the type of higher order moment:
%           (i) 'skew1', skewness
%           (ii) 'skew2', skewness correcting for bias
%           (iii) 'kurt1', kurtosis
%           (iv) 'kurt2', kurtosis correcting for bias
% 

function out = DN_Cumulants(y,whatcum)
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
    error('Unknown cumulant ''%s'' specified: should be ''skew1'', ''skew2'', ''kurt1'', or ''kur2''',whatcum)
end

end