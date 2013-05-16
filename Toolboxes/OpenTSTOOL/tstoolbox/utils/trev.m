function result=trev(x,tau)

% TREV calculates the time reversibility coefficient
% belongs to the surrogate data test statistics.
%
% x - column vector
% tau - delay in samples


%x1=data(cut(s,1,1+tau,length(data(s))))';
x1=x(1:end-tau);
if length(x(:,1))~=1
  error('TREV can only manage scalar time series, sorry...');
end

%x2=data(cut(s,1,1,length(data(s))-tau))';
x2=x(1+tau:end);
up=(x1-x2).^3;
down=(x1-x2).^2;

result=mean(up)/(mean(down)^(3/2));

return
 