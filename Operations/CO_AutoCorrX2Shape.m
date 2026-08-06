function out = CO_AutoCorrX2Shape(y, maxLag)
% CO_AutoCorrX2Shape   Shape of the time-reversibility profile of a time series.
%
% CO_AutoCorrX2 computes two asymmetric, 'leverage'-type lag-profiles:
%
%   forward(tau):  <x_t.x_{t+tau}^2>   (signed value now, energy later)
%   backward(tau): <x_t^2.x_{t+tau}>   (energy now, signed value later)
%
% For a time-reversible process these coincide at every lag (any shared linear
% correlation structure contributes equally to both); a systematic difference,
%   diff(tau) = forward(tau) - backward(tau),
% is therefore a lag-resolved time-irreversibility statistic, generalizing the
% single-lag trev/tc3-style statistics to a full profile, cf. the
% leverage-effect correlation function of Bouchaud, Matacz & Potters, Phys.
% Rev. Lett. 87, 228701 (2001).
%
% This function characterizes the SHAPE of diff(tau) across lags -- its decay,
% persistence, and extrema -- mirroring how CO_AutoCorrShape characterizes the
% shape of the ordinary ACF. (An earlier version of this function instead
% characterized the forward and backward profiles' shapes separately, but on
% 300 real time series from INP_Empirical1000.mat their shape descriptors
% (centroid decay timescale, profile-autocorrelation, etc.) were correlated at
% r=0.84-0.97 with each other -- i.e., overwhelmingly redundant, since both
% profiles inherit most of their shape from whatever ordinary linear
% correlation the series has. The difference profile cancels that shared
% component and isolates the genuinely asymmetric/nonlinear structure.)
%
% ---INPUTS:
% y, the input time series (should be z-scored: zero mean, unit variance)
% maxLag, the maximum lag to compute the profile up to. Can be a positive
%       integer, or the string 'doubleDrown', which sets it to twice the
%       first zero-crossing of the ordinary (linear) autocorrelation
%       function (cf. the 'doubleDrown' option of CO_AutoCorrShape), bounded
%       to lie in [10, floor(N/4)].

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Check inputs & set defaults:
% ------------------------------------------------------------------------------
y = y(:);
N = length(y);

if nargin < 2 || isempty(maxLag)
	maxLag = 'doubleDrown';
end

if ischar(maxLag)
	switch maxLag
		case 'doubleDrown'
			tau0 = CO_FirstCrossing(y, 'ac', 0, 'discrete');
			if isnan(tau0) || tau0 == 0
				tau0 = 10; % fallback for pathological/near-constant series
			end
			maxLag = min(floor(N / 4), max(10, 2 * tau0));
		otherwise
			error('Unknown maxLag setting: ''%s''', maxLag);
	end
end

if maxLag < 5
	% Too short a series/profile to say anything meaningful
	out = NaN;
	return
end

% ------------------------------------------------------------------------------
%% Compute the forward and backward lag-profiles, and their difference:
% ------------------------------------------------------------------------------
gFwd = CO_AutoCorrX2(y, 1:maxLag, 'forward');
gBwd = CO_AutoCorrX2(y, 1:maxLag, 'backward');

if any(isnan(gFwd)) || any(isnan(gBwd))
	out = NaN;
	return
end

diffProfile = gFwd(:) - gBwd(:);

% -------------------------------------------------------------------------------
% Lag-1 difference: the single-lag time-irreversibility statistic, comparable
% to CO_trev. (Note: the raw lag-1 and lag-0 values themselves -- gFwd(1),
% gBwd(1), and mean(y.^3) -- are not included as outputs here since they
% duplicate existing operations exactly: AC_nl_01, AC_nl_11, and DN_Moments_3
% respectively.)
% -------------------------------------------------------------------------------
out.diff1 = diffProfile(1);

% -------------------------------------------------------------------------------
% Basic stats on the difference profile
% -------------------------------------------------------------------------------
out.sumdiff = sum(diffProfile);
out.meandiff = mean(diffProfile);
out.meanabsdiff = mean(abs(diffProfile));
out.rmsdiff = sqrt(mean(diffProfile.^2));

% -------------------------------------------------------------------------------
% Characteristic (centroid) decay timescale of the difference profile -- how
% many lags the time-irreversibility signature persists over. Centroid-based
% (rather than an exponential fit) to avoid a curve-fitting-toolbox
% dependency and to remain well-defined for non-monotonic profiles.
% -------------------------------------------------------------------------------
tau = (1:maxLag)';
if sum(abs(diffProfile)) > 0
	out.centroiddiff = sum(tau .* abs(diffProfile)) / sum(abs(diffProfile));
else
	out.centroiddiff = NaN;
end

% -------------------------------------------------------------------------------
% Autocorrelation of the difference profile (smoothness/persistence of the
% irreversibility signature itself), cf. the ac1 field of CO_AutoCorrShape
% -------------------------------------------------------------------------------
out.ac1diff = CO_AutoCorr(diffProfile, 1, 'Fourier');

% -------------------------------------------------------------------------------
% Local extrema of the difference profile, cf. CO_AutoCorrShape
% -------------------------------------------------------------------------------
ddiff = diff(diffProfile);
dddiff = diff(ddiff);
extrr = BF_SignChange(ddiff, 1);
sdsp = dddiff(extrr);
out.nminima = sum(sdsp > 0);
out.nmaxima = sum(sdsp < 0);
out.pextrema = length(sdsp) / maxLag;

% -------------------------------------------------------------------------------
% How quickly the time-irreversibility signature changes sign
% -------------------------------------------------------------------------------
signChangeIdx = find(sign(diffProfile(1:end - 1)) ~= sign(diffProfile(2:end)), 1, 'first');
if isempty(signChangeIdx)
	out.firstsignchangediff = maxLag; % no sign change within the window measured
else
	out.firstsignchangediff = signChangeIdx;
end

% -------------------------------------------------------------------------------
% Shape similarity of the two profiles (a different angle from their
% difference: do they have a similar shape regardless of overall magnitude?)
% -------------------------------------------------------------------------------
cc = corrcoef(gFwd(:), gBwd(:));
out.corrfwdbwd = cc(1, 2);

out.maxLag = maxLag; % record how far the profile was measured

end
