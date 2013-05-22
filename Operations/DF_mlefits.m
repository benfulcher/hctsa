function out = DF_mlefits(y,fitwhat,outwhat)
% Uses the mle function from Matlab's Statistics Toolbox to fit distributions to the time series, y
% INPUTS:
% fitwhat specifies which fit to do
% outwhat specifies which component of the output to take
% Ben Fulcher 2009

switch fitwhat
	case 1 % Gaussian Fit
		phat = mle(y);
		out = phat(outwhat);
		% outwhat = 1: mean of gaussian fit
		% outwhat = 2: std of gaussian fit
	case 2
		phat = mle(y,'distribution','uniform');
		out = phat(outwhat); % outwhat = (1 = a), (2 = b)
	case 3
		phat = mle(y,'distribution','geometric');
		out = phat(outwhat);
    otherwise
        error('DF_mlefits: Invalid fit specified')
end

end