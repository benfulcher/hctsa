function out = CO_fzcglscf(y,alpha,beta)
% Calculates the generalized self-correlation function, as introduced by
% Duarte Queiros and Moyano in Physica A, 2007
% Keeps calculating until crosses zero, returns this lag (maximum = 400)

glscfs=zeros(400,1);

for i=1:400
	tau = i;
	y1 = abs(y(1:end-tau));
	y2 = abs(y(1+tau:end));

	glscfs(i) = (mean((y1.^alpha).*(y2.^beta)) - mean(y1.^alpha)*mean(y2.^beta)) / ...
	 		( sqrt(mean(y1.^(2*alpha)) - mean(y1.^alpha)^2) * sqrt(mean(y2.^(2*beta)) - mean(y2.^beta)^2) );

	if i>1 && glscfs(i)*glscfs(i-1)<0,
		% draw a straight line between these two and look at where hits zero
		out = i-1 + glscfs(i)/(glscfs(i)-glscfs(i-1));
		return;
	end
end

out=i;


end