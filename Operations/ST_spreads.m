function out = ST_spreads(y,meth)
% time series 'y'
% spread method 'meth'
% Ben Fulcher, 2008

if nargin < 2 || isempty(meth)
    meth = 'std'; % return std by default
end

switch meth
	case 'std' % standard deviation
		out = std(y);
	case 'iqr' % interquartile range
		out = iqr(y);
	case 'mad' % mean absolute deviation
		out = mad(y,0);
    case 'mead' % median absolute deviation
        out = mad(y,1);
    otherwise
        error('Unknown spread measure ''%s''',meth)
end

end