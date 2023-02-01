function f = CR_RAD(x,tau,doAbs)
% RAD Compute the rescaled auto-density, a metric for inferring the
% distance to criticality that is insensitive to uncertainty in the noise strength.
% Calibrated to experiments on the Hopf bifurcation with variable and unknown
% measurement noise.
%
% Devised and implemented by Brendan Harris, @brendanjohnharris (GitHub), 2023.

%---INPUTS:
%   x:        The input time series (vector).
%   doAbs:    Whether to centre the time series at 0 then take absolute values (logical flag)
%   tau:      The embedding and differencing delay in units of the timestep (integer)
%
%---OUTPUTS:
%   f:        The RAD feature value

%-------------------------------------------------------------------------------
% Check inputs, set defaults
%-------------------------------------------------------------------------------
if nargin < 2 || isempty(tau)
    tau = 1;
end
if nargin < 3 || isempty(doAbs)
    doAbs = true;
end
%-------------------------------------------------------------------------------
% Basic checks & preprocessing
%-------------------------------------------------------------------------------
if isrow(x)
    x = x';
end
if doAbs
    x = x - median(x);
    x = abs(x);
end
if ischar(tau) && strcmp(tau,'tau')
    % Make tau the first zero crossing of the autocorrelation function
    tau = CO_FirstCrossing(x,'ac',0,'discrete');
end
%-------------------------------------------------------------------------------

% Delay embed at interval tau, m = 2
y = x(tau+1:end);
x = x(1:end-tau);

% Median split
subMedians = (x < median(x));
superMedianSD = std(x(~subMedians));
subMedianSD = std(x(subMedians));

% Properties of the auto-density
sigma_dx = std(y - x);
densityDifference = 1./superMedianSD - 1./subMedianSD;

f = sigma_dx.*densityDifference;

end
