function g = gauss(N)

% gauss(N)
% 
% returns N normally distributed random numbers
% with zeo mean and 
% reference : Num. Recipes, Chapter 7.2 Normal Deviates

rsq  = [];

M = ceil(N/2);
while length(rsq) < M	% make shure we really have at least M values
	v = 2*rand(ceil(M*1.33),2)-1;		% produce more random numbers
	rsq = v(:,1).*v(:,1)+v(:,2).*v(:,2);
	ind = find((rsq >=1) | (rsq == 0));
	rsq(ind) = [];						% because we want to remove some
end

v(ind,:) = [];
v = v(1:M,:);
rsq = rsq(1:M);
fac = sqrt(-2 * log(rsq) ./ rsq);
g = [v(:,1) .* fac ; v(:,2) .* fac];
g = g(1:N);
