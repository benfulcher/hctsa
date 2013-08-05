% CO_acfdiff
% 
% Estimates the autocorrelation function using a fast Fourier Transform method
% implemented in TSTOOL and returns the mean square discrepancy between the
% autocorrelation coefficients obtained in this way from those obtained in the
% time domain using CO_AutoCorr.
% 
% TSTOOL: http://www.physik3.gwdg.de/tstool/
% 
% No real rationale behind this, other than the difference in autocorrelations
% computed by the two methods may somehow be informative of something about the
% time series...
% 

function out = CO_acfdiff(y,maxlag)
% Ben Fulcher, October 2009

if nargin < 2 || isempty(maxlag)
    maxlag = 50; % compare across the first maxlag autocorrelations
end

% First maxlag autocorrelations
co_fft = data(acf(signal(y),maxlag*2));

nlags = length(co_fft);
co_ben = zeros(nlags,1);
for i = 1:nlags
    co_ben(i) = CO_AutoCorr(y,i-1);
end

out = norm(co_ben - co_fft);

end