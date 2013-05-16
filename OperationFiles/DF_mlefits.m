function out = DF_mlefits(y,n,i)
% n specifies which fit to do
% i specifies which component of the output to take
% Ben Fulcher 2009

switch n
	case 1 % Gaussian Fit
		phat = mle(y);
		out = phat(i);
		% i=1: mean of gaussian fit
		% i=2: std of gaussian fit
	case 2
		phat = mle(y,'distribution','uniform');
		out = phat(i); % i = (1 = a), (2 = b)
	case 3
		phat = mle(y,'distribution','geometric');
		out = phat(i);
end
