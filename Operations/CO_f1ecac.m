function out = CO_f1ecac(y)
% Calculates the 1/e correlation length of the time series using RM_information
% Ben Fulcher 2008
  
oone = 1/exp(1);
a(1) = RM_information(y,y); % very weird -- why is this not 1?? Or use autocor?

for i = 2:length(y)-1
    a(i) = CO_autocorr(y,i);
    if (a(i-1)-oone)*(a(i)-oone) < 0
        out = i;
        return
    end
end

% If no minimum in entire spectrum return the maximum value
out = length(y);

end