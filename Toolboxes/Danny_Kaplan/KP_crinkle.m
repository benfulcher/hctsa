function out = KP_crinkle(x)
% Calculates the "crinkle statistic" on a vector x
% 	<(x_{t-1}-2*x_t+x_{t+1})^4> / < ( x_t^2 ) >^2
%	as proposed by James Theiler
% Copyright (c) 1996 by D. Kaplan, All Rights Reserved

% Subtract out the mean 
x = x - mean(x);
x2 = mean(x.*x)^2;

if x2 == 0 
	out = 0;
	return
end

d2 = 2*x(2:end-1) - x(1:end-2) - x(3:end);
out = mean(d2.^ 4)/x2;

end
