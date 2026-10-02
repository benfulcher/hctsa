function f = CR_RAD(x, tau, doAbs)
% CR_RAD   The rescaled auto-density, a metric for inferring the distance to criticality.
%
% Computes the rescaled auto-density (RAD), a metric for inferring the distance to
% criticality that is insensitive to uncertainty in the noise strength. Calibrated to
% experiments on the Hopf bifurcation with variable and unknown measurement noise.
%
% With doAbs, the series is first centered at its median and made non-negative. It is
% then delay-embedded in two dimensions, (x(t), x(t+tau)), and split at the median of
% x(t). The output is the standard deviation of the increments x(t+tau) - x(t), times
% (1/std of the upper half of x(t) - 1/std of the lower half of x(t)).
%
% ---INPUTS:
% x, the input time series (vector)
% tau, the embedding and differencing delay in units of the timestep (integer; default
%      1; 'tau' sets it to the first zero-crossing of the autocorrelation function of
%      the series after any doAbs transformation)
% doAbs, whether to center the time series at its median and then take absolute values
%        (logical flag; default true)
%
% ---OUTPUTS:
% f, a scalar: the RAD feature value.
%
% ---REFERENCES:
% B. Harris, L. L. Gollo and B. D. Fulcher, "Tracking the distance to criticality in
% systems with unknown noise", Physical Review X 14(3), 031021 (2024).
% DOI: 10.1103/PhysRevX.14.031021
% Please cite this paper if you use this function in your work.
%
% ---NOTES:
% Devised and authored by Brendan Harris, 2023 (@brendanjohnharris on GitHub).
% Edits by Ben Fulcher for incorporating into hctsa.

% -------------------------------------------------------------------------------
% Check inputs, set defaults
% -------------------------------------------------------------------------------
if nargin < 2 || isempty(tau)
	tau = 1;
end
if nargin < 3 || isempty(doAbs)
	doAbs = true;
end
% -------------------------------------------------------------------------------
% Basic checks & preprocessing
% -------------------------------------------------------------------------------
if isrow(x)
	x = x';
end
if doAbs
	x = x - median(x);
	x = abs(x);
end
if ischar(tau) && ismember(tau, {'ac1e', 'mi'})
	% Adaptive delay: see BF_GetTau
	tau = BF_GetTau(x, tau);
	if isnan(tau)
		f = NaN; return
	end
end
if ischar(tau) && strcmp(tau, 'tau')
	% Make tau the first zero crossing of the autocorrelation function
	tau = CO_FirstCrossing(x, 'ac', 0, 'discrete');
	if isnan(tau)
		f = NaN; return
	end
end
% -------------------------------------------------------------------------------

% Delay embed at interval tau, m = 2
y = x(tau + 1:end);
x = x(1:end - tau);

% Median split
subMedians = (x < median(x));
superMedianSD = std(x(~subMedians));
subMedianSD = std(x(subMedians));

% Properties of the auto-density
sigma_dx = std(y - x);
densityDifference = 1 ./ superMedianSD - 1 ./ subMedianSD;

f = sigma_dx .* densityDifference;

end
