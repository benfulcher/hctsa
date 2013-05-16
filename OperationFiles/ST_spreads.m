function out = ST_spreads(y,meth)
% time series 'y'
% spread method 'meth'

switch meth
	case 'std' % standard deviation
		out = std(y);
	case 'iqr' % interquartile range
		out = iqr(y);
	case 'mad' % mean absolute deviation
		out = mad(y,0);
    case 'mead' % median absolute deviation
        out = mad(y,1);
%     case 'RMSSD'
%         % a simple statistical HRV measure
    otherwise
        error('Invalid spread measure specified')
end

end