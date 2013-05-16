
function out=MF_ARfits(y,p,n)
% AR fit of order p
% obtain 'nth' statistic

[a e]=arcov(y,p);

if n==0
	out=e; % variance of fit...()?)
else
	out=a(n);
end
end