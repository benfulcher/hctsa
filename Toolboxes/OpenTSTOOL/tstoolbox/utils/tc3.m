function result=tc3(x,tau)

% TC3 calculates the higher moment coefficient
% belongs to the surrogate data test statistics.
%
% x - column vector
% tau - delay in samples

%x=data(s);

if length(x(:,1))~=1
  error('TC3 can only manage scalar time series, sorry...');
end

up=x(1+2*tau:end).*x(1+tau:end-tau).*x(1:end-2*tau);
down=x(1+2*tau:end).*x(1+tau:end-tau);


result=mean(up)/(abs(mean(down))^(3/2));


return
