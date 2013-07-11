function out = CO_fzcmi(y)
% Find the first minimum of the automutual information calculated using RM_information for an input time series, y
% Ben Fulcher 19/3/2010 -- corrected for error and cleaned up a few things
% from this old, messy code (probably from 2008). In fact, I think it's used only as CO_fmmi?

N = length(y); % the length of the time series

a = zeros(N-3,1);
a(1) = RM_information(y,y);
for i = 2:N-3
    try
        a(i) = RM_information(y(1:end-i),y(1+i:end));
    catch emsg
        % Some crazy thing happened
        disp('Can you believe that some crazy thing happened with CO_fzcmi??')
        disp(['Well it did: ' emsg]);
        out = NaN;
        return
    end
    if a(i-1)*a(i) < 0
        out = i;
        return
    end
end

% If no minimum in entire spectrum
out = N;

end