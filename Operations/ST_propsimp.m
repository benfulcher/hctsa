function out = ST_propsimp(x,condition)
% Returns the proportion of the time series that is a given value
% i.e., zeros, positive, and values that are greater than or equal to zero
% Ben Fulcher, 2009

N = length(x); % length of the time series

switch condition
    case 'zeros' % returns the proportion of zeros in the input vector
        out = sum(x == 0)/N;
    case 'positive'
        out = sum(x > 0)/N;
    case 'geq0'
        out = sum(x >= 0)/N;
    otherwise
        error('Unknown condition ''%s''',condition);
end

end