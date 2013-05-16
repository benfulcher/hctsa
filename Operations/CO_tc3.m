function out = CO_tc3(y,tau)
% Implements the TC3 function (which I got from the TSTOOL package)
% Ben Fulcher 15/11/2009


if strcmp(tau,'ac')
    tau = CO_fzcac(y);
elseif strcmp(tau,'mi')
    tau = CO_fmmi(y);
end

yn = y(1:end-2*tau);
yn1 = y(1+tau:end-tau); % yn1, tau steps ahead
yn2 = y(1+2*tau:end); % yn2, 2*tau steps ahead


% The expression used in TSTOOL tc3
out.raw = mean(yn.*yn1.*yn2)/abs(mean(yn.*yn1))^(3/2);

% The magnitude
out.abs = abs(out.raw);

% The numerator
out.num = mean(yn.*yn1.*yn2);
out.absnum = abs(out.num);

% The denominator
out.denom = abs(mean(yn.*yn1))^(3/2);

end