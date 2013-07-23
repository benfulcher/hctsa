% DN_ProportionValues
% 
% Returns statistics on the values of the raw time series: the proportion
% of zeros in the raw time series, the proportion of positive values, and the
% proportion of values greater than or equal to zero.
% 
% INPUTS:
% x, the input time series
% 
% propwhat, the proportion of a given type of value in the time series:
%           (i) 'zeros': values that equal zero
%           (ii) 'positive': values that are strictly positive
%           (iii) 'geq0': values that are greater than or equal to zero
% 

function out = DN_ProportionValues(x,propwhat)
% Ben Fulcher, 2009

N = length(x); % length of the time series

switch propwhat
    case 'zeros' % returns the proportion of zeros in the input vector
        out = sum(x == 0)/N;
        
    case 'positive'
        out = sum(x > 0)/N;
        
    case 'geq0'
        out = sum(x >= 0)/N;
        
    otherwise
        error('Unknown condition to measure: ''%s''',propwhat);
end

end