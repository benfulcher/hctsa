% DF_mlefits
% 
% Fits either a Gaussian, Uniform, or Geometric distribution to the data using
% maximum likelihood estimation via the MATLAB function mle
% from the Statistics Toolbox.
% 
% INPUTS:
% y, the time series
% fitwhat, the type of fit to do: 'gaussian', 'uniform', or 'geometric'.
% 
% Outputs are parameters from the fit.

function out = DF_mlefits(y,fitwhat)
% Ben Fulcher, 2009

if nargin < 2
    fitwhat = 'gaussian'; % fit a Gaussian by default
end

switch fitwhat
case 'gaussian'
	phat = mle(y);
	out.mean = phat(1); % mean of Gaussian fit
    out.std = phat(2); % std of Gaussian fit
case 'uniform' % turns out to be shit
    phat = mle(y,'distribution','uniform');
    out.a = phat(1);
    out.b = phat(2);
case 'geometric'
    out = mle(y,'distribution','geometric'); % just a single output
otherwise
    error('Invalid fit specifier, ''%s''',fitwhat)
end

end