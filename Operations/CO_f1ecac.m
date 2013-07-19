% CO_f1ecac
% 
% Finds where autocorrelation function first crosses 1/e, the 1/e correlation
% length
% 
% INPUTS:
% y, the input time series

function out = CO_f1ecac(y)
% Ben Fulcher, 2008
  
N = length(y); % time-series length
oone = 1/exp(1);

for i = 1:N-1
    a(i) = CO_autocorr(y,i);
    if (i > 1) && ((a(i-1)-oone)*(a(i)-oone) < 0)
        % Crossed the 1/e line
        out = i;
        return
    end
end

% If no minimum in entire spectrum return the maximum value
out = N;

end