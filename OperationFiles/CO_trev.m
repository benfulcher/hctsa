function out = CO_trev(y,tau)
% Implements the trev function (which I got from the TSTOOL package, and elsewhere...)
% Ben Fulcher 15/11/2009


if strcmp(tau,'ac')
    tau = CO_fzcac(y);
elseif strcmp(tau,'mi')
    tau = CO_fmmi(y);
end

yn = y(1:end-tau);
yn1 = y(1+tau:end); % yn, tau steps ahead

% keyboard

% The expression used in TSTOOL
out.raw = mean((yn1-yn).^3)/(mean((yn1-yn).^2))^(3/2);

% The magnitude
out.abs = abs(out.raw);

% The numerator
out.num = mean((yn1-yn).^3);
out.absnum = abs(out.num);

% The denominator
out.denom = (mean((yn1-yn).^2))^(3/2);

end