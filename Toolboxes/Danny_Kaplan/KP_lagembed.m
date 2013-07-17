function y = KP_lagembed(x,M,lag)
% lagEmbed(x,dim,lag) constructs an embedding of a time series on a vector
% lagEmbed(x,dim) makes an m-dimensional embedding with lag 1
% lagEmbed(x,dim,lag) uses the specified lag
% Copyright (c) 1996 by D. Kaplan, All Rights Reserved

if nargin < 3
	lag = 1;
end
if nargin < 4
	advance=0;
end

%convert x to a column
[xr,xc] = size(x);
if xr == 1	
    x = x';
end


lx = length(x);
	
newsize = lx - lag*(M-1);
y = zeros(newsize,M);
i=1;

for j = 0:-lag:lag*(-(M-1))

	first=1+lag*(M-1)+j;

	last=first+newsize-1;


	if last > lx

		last = lx;

	end

	y(:,i) = x(first:last, 1);

	i = i+1;

end

