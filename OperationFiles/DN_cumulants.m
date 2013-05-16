function out=DN_cumulants(y,n)
switch n
	case 'skew1' % skewness
		out=skewness(y);
    case 'skew2' % corrects for bias
        out=skewness(y,0);
	case 'kurt1' % kurtosis
		out=kurtosis(y);
    case 'kurt2' % corrects for bias
        out=kurtosis(y,0);
end