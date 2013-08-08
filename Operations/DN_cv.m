% DN_cv
% 
% Calculates the coefficient of variation, sigma^k / mu^k, of order k
% 
% INPUTS:
% 
% x, the input time series
% 
% k, the order of coefficient of variation (k = 1 is usual)
% 

function out = DN_cv(x,k)
% Ben Fulcher, 2009

% Check inputs
if nargin < 2 || isempty(k)
    k = 1; % Do standard CV by default
end
if (rem(k,1) ~= 0) || (k < 0)
    warning('k should probably be a positive integer');
    % Carry on with just this warning, though
end

out = (std(x))^k / (mean(x))^k;

end