% DN_Spread
% 
% Returns the spread of the raw time series, as the standard deviation,
% inter-quartile range, mean absolute deviation, or median absolute deviation.
% 
% INPUTS:
% y, the input time series
% 
% thespread, the spead measure:
%               (i) 'std': standard deviation
%               (ii) 'iqr': interquartile range
%               (iii) 'mad': mean absolute deviation
%               (iv) 'mead': median absolute deviation
% 

function out = DN_Spread(y,thespread)
% Ben Fulcher, 2008

if nargin < 2 || isempty(thespread)
    thespread = 'std'; % return std by default
end

switch thespread
	case 'std' % standard deviation
		out = std(y);
        
	case 'iqr' % interquartile range
		out = iqr(y);
        
	case 'mad' % mean absolute deviation
		out = mad(y,0);
        
    case 'mead' % median absolute deviation;
        out = mad(y,1);
        
    otherwise
        error('Unknown spread measure ''%s''',thespread)
end

end