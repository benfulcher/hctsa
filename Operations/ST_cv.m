function out = ST_cv(x,k)
% Calculates the coefficient of variation, CV
% k a positive integer: the 'moment' of cv
% Ben Fulcher, 2009

% Check inputs
if nargin < 2
    k = 1;
    % Do standard CV by default
end
if rem(k,1) ~= 0 || k < 0
    warning('k should be a positive integer');
    % Carry on with just this warning, though
end

out = (std(x))^k / (mean(x))^k;

end