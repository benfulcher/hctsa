function out = DF_mlefits(y,fitwhat)
% Uses the mle function from Matlab's Statistics Toolbox to fit distributions to the time series, y
% INPUTS:
% fitwhat specifies which fit to do
% outwhat specifies which component of the output to take
% Ben Fulcher 2009
% Ben Fulcher 2013 -- turned into a master operation with sensible inputs

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

% OLD VERSION:
% switch fitwhat
%     case 1 % Gaussian Fit
%         outwhat = 1: mean of gaussian fit
%         outwhat = 2: std of gaussian fit
%     case 2
%         phat = mle(y,'distribution','uniform');
%         out = phat(outwhat); % outwhat = (1 = a), (2 = b)
%     case 3
%         phat = mle(y,'distribution','geometric');
%         out = phat(outwhat);
%     otherwise
%         error('DF_mlefits: Invalid fit specified')
% end

end