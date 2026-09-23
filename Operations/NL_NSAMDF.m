function out = NL_NSAMDF(x, tauMult, winLenRel, shiftLenRel, degree, doPlot)
% NL_NSAMDF computes the nonlinearity measure L through nsAMDF
% (nonlinear average magnitude difference function), developed by:
%
% Ozkurt et al. (2020), "Identification of nonlinear features in cortical and
% subcortical signals of Parkinson's Disease patients via a novel efficient
% measure", NeuroImage.
%
% Please refer to and cite this paper if you use this function in your work.
%
% The nsAMDF of degree p is the p-norm of the lag-tau increments,
% ||x(t+tau) - x(t)||_p, as a function of tau = 0, ..., maxLag, averaged over
% sliding windows. For a linear Gaussian process the increments at every lag are
% Gaussian, so the normalized p = 2 and p = degree curves have the same shape;
% L measures how far they differ, which is sensitive to non-sinusoidal waveform
% shape and other departures from Gaussian increment structure. It does not
% detect all nonlinearity (e.g., a threshold AR with near-Gaussian increments).
%
% The paper used a lag range of 1 s and a window of 14 s at the sampling rate of
% the data. hctsa has no sampling rate, so both are set relative to the
% correlation timescale of the series (first zero-crossing of the
% autocorrelation function, tau_ac): maxLag = ceil(tauMult*tau_ac), and the
% window length is winLenRel*maxLag.
%
% ---INPUTS:
%
% x, the input time series.
%
% tauMult, the maximum lag as a multiple of tau_ac.
%
% winLenRel, the window length as a multiple of the maximum lag.
%
% shiftLenRel, the window shift as a proportion of the window length.
%
% degree, the chosen degree p (> 2) should ideally be large enough to capture the
%           highest order of nonlinearity within the data.
%           The paper used p = 7 for Parkinsonian data.
%
% doPlot, true to plot nsAMDF sequences.
%
% ---OUTPUTS:
%
% L: root-mean-square difference between the normalized nsAMDF curves of degree
%       2 and of degree p, across lags 0, ..., maxLag (scale-invariant).
%
% s2: normalized nsAMDF for the degree 2
%
% sd: normalized nsAMDF for the chosen degree greater than 2
%
% Required subfunctions are NormedSingleCurveLengthWindowed.m & NormedSingleCurveLength.m
%
%   Authored by Tolga Esat Ozkurt, 2020. (tolgaozkurt@gmail.com)
%   Edits by Ben Fulcher for incorporating into hctsa.

% -------------------------------------------------------------------------------
% Set defaults:
if nargin < 2 || isempty(tauMult)
	tauMult = 2;
end
if nargin < 3 || isempty(winLenRel)
	winLenRel = 10;
end
if nargin < 4 || isempty(shiftLenRel)
	shiftLenRel = 0.5;
end
if nargin < 5 || isempty(degree)
	degree = 7;
end
if nargin < 6
	doPlot = false;
end

% -------------------------------------------------------------------------------
% Set the lag range and window from the correlation timescale:
tau = CO_FirstCrossing(x, 'ac', 0, 'discrete');
if isnan(tau)
	out = NaN; return % data-dependent: no correlation length could be estimated
end
maxLag = ceil(tauMult * tau);
windowLength = winLenRel * maxLag;
if windowLength > length(x)
	out = NaN; return % data-dependent: too short relative to its correlation time
end
shiftLength = max(1, floor(shiftLenRel * windowLength));

% -------------------------------------------------------------------------------
% nsAMDF for p = 2:
s2 = NormedSingleCurveLengthWindowed(x, windowLength, shiftLength, maxLag, 1, 2);
out.s2 = s2 ./ max(s2); % normalized

% nsAMDF for p = degree:
sd = NormedSingleCurveLengthWindowed(x, windowLength, shiftLength, maxLag, 1, degree);
out.sd = sd ./ max(sd); % normalized

% Compare the normalized curves (the original code compared the unnormalized
% ones, which makes L scale with the amplitude of the data):
out.L = sqrt(mean((out.s2 - out.sd).^2));

% -------------------------------------------------------------------------------
if doPlot
	figure
	plot(0:maxLag, out.s2, 'b')
	hold on
	plot(0:maxLag, out.sd, 'g')
	xlabel('Lag'); legend('p = 2', sprintf('p = %u', degree))
end

end
