% DN_MinMax
% 
% Returns the maximum and minimum values of the input time series.
% 
% INPUTS:
% 
% y, the input time series
% 
% minormax, either 'min' or 'max' to return either the minimum or maximum of y
% 

function out = DN_MinMax(y,minormax)
% Ben Fulcher, 2008

switch minormax
    case 'max'
        out = max(y);
        
    case 'min'
        out = min(y);
        
    otherwise
        error('Unknown method ''%s''',minormax)
end

end