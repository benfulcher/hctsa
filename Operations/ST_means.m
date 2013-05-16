function out = ST_means(y,meantype)
N = length(y);

switch meantype
	case 'norm' % mean
		out = mean(y);
	case 'geom' % geometric mean
		out = geomean(y); %(prod(y))^(1/N);
	case 'harm' % harmonic mean
		out = harmmean(y); %N/sum(y.^(-1));
	case 'rms' % rms
		out = sqrt(sum(y.^2)/N);
    case 'iqm' % interquartile mean
        p = prctile(y, [25; 75]);
        out = mean(y(y>p(1) & y<p(2)));
    case 'midhinge' % average of 1st and third quartiles
        p = prctile(y, [25; 75]);
        out = mean(p);
end

end