% DN_Mean
% 
% Measures a given type of 'mean', or measure of location of the time series.
% 
% INPUTS:
% 
% y, the input time series
% 
% meantype, (i) 'norm' or 'arithmetic', arithmetic mean
%           (ii) 'median', median
%           (iii) 'geom', geometric mean
%           (iv) 'harm', harmonic mean
%           (v) 'rms', root-mean-square
%           (vi) 'iqm', interquartile mean
%           (vii) 'midhinge', midhinge
% 

function out = DN_Mean(y,meantype)
% Ben Fulcher, 2008

% Check Inputs
if nargin < 2 || isempty(meantype)
    meantype = 'arithmetic'; % normal mean
end

N = length(y); % time-series length

switch meantype
	case {'norm','arithmetic'} % mean
		out = mean(y);
        
    case 'median' % median
        out = median(y);
        
	case 'geom' % geometric mean
		out = geomean(y); %(prod(y))^(1/N);
        
	case 'harm' % harmonic mean
		out = harmmean(y); %N/sum(y.^(-1));
        
	case 'rms' % rms
		out = sqrt(sum(y.^2)/N);
        
    case 'iqm' % interquartile mean
        p = prctile(y, [25; 75]);
        out = mean(y(y >= p(1) & y <= p(2)));
        
    case 'midhinge' % average of 1st and third quartiles
        p = prctile(y, [25; 75]);
        out = mean(p);
        
    otherwise
        error('Unknown mean type ''%s''', meantype);
end

end