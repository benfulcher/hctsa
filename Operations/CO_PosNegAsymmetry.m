function out = CO_PosNegAsymmetry(y)
% CO_PosNegAsymmetry   Asymmetry of local dynamics between positive and negative regimes.
%
% Splits the time series by the sign of each value (assumes y is z-scored, so
% the split is around the mean) and asks whether the one-step-ahead dynamics
% differ between the two regimes: is the series more volatile, or more
% persistent (higher one-step autocorrelation), when the current value is
% above vs. below its mean? This targets a form of distribution-dynamics
% interaction not captured by CO_StickAngles, which instead compares the
% distribution of local slopes *within* each same-sign subsequence.
%
% ---INPUTS:
% y, the input time series (assumed z-scored: the regime split threshold is 0)
%
% ---OUTPUTS: the conditional volatility and one-step autocorrelation of the
%            positive/negative regimes, and normalized contrasts between
%            them (the volatility contrast is a leverage-effect-style
%            statistic; the autocorrelation contrast is a threshold-AR(1)-
%            style statistic). Also isolates the two zero-crossing
%            transition types (positive-to-negative, negative-to-positive)
%            -- posMask/negMask are conditioned on the current value only,
%            so they mix crossing and non-crossing steps together; the
%            crossing-specific fields ask instead whether the jump *at* a
%            regime switch is itself asymmetric (e.g. sharper downward
%            crossings than upward ones).

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

y = y(:);

% ------------------------------------------------------------------------------
% Regime-conditioned one-step-ahead pairs: yLag(t) = y(t), yNext(t) = y(t+1)
% ------------------------------------------------------------------------------
yLag = y(1:end - 1);
yNext = y(2:end);
dy = diff(y);

posMask = (yLag >= 0);
negMask = ~posMask;
nPos = sum(posMask);
nNeg = sum(negMask);

out.propPos = mean(posMask); % proportion of time spent in the positive regime

% ------------------------------------------------------------------------------
% Conditional volatility (leverage-effect-style asymmetry)
% ------------------------------------------------------------------------------
if nPos >= 2
	out.volPos = std(dy(posMask));
else
	out.volPos = NaN;
end
if nNeg >= 2
	out.volNeg = std(dy(negMask));
else
	out.volNeg = NaN;
end
volAll = std(dy);
if nPos >= 2 && nNeg >= 2 && volAll > 0
	% Normalized contrast: positive when the positive regime is more volatile.
	out.volAsym = (out.volPos - out.volNeg) / volAll;
else
	out.volAsym = NaN;
end

% ------------------------------------------------------------------------------
% Conditional one-step autocorrelation (threshold-AR(1)-style asymmetry)
% ------------------------------------------------------------------------------
if nPos >= 2
	R = corrcoef(yLag(posMask), yNext(posMask));
	out.ac1Pos = R(2, 1);
else
	out.ac1Pos = NaN;
end
if nNeg >= 2
	R = corrcoef(yLag(negMask), yNext(negMask));
	out.ac1Neg = R(2, 1);
else
	out.ac1Neg = NaN;
end
% Difference on the (already bounded, comparable) correlation scale:
out.ac1Asym = out.ac1Pos - out.ac1Neg;

% ------------------------------------------------------------------------------
% Zero-crossing transitions: positive-to-negative (PN) and negative-to-
% positive (NP). Unlike posMask/negMask (conditioned on the current value
% only), these isolate the step *at* which the regime actually switches.
% ------------------------------------------------------------------------------
pnMask = posMask & (yNext < 0); % downward crossing
npMask = negMask & (yNext >= 0); % upward crossing
nPN = sum(pnMask);
nNP = sum(npMask);

out.propPN = mean(pnMask); % proportion of steps that are downward crossings
out.propNP = mean(npMask); % proportion of steps that are upward crossings

% Crossing-jump volatility:
if nPN >= 2
	out.volPN = std(dy(pnMask));
else
	out.volPN = NaN;
end
if nNP >= 2
	out.volNP = std(dy(npMask));
else
	out.volNP = NaN;
end
if nPN >= 2 && nNP >= 2 && volAll > 0
	% Positive when downward crossings are more violent than upward ones.
	out.volAsymCross = (out.volPN - out.volNP) / volAll;
else
	out.volAsymCross = NaN;
end

% Does the pre-crossing value predict the post-crossing value (i.e. does a
% deeper excursion before crossing predict a deeper overshoot after it)?
if nPN >= 2
	R = corrcoef(yLag(pnMask), yNext(pnMask));
	out.ac1PN = R(2, 1);
else
	out.ac1PN = NaN;
end
if nNP >= 2
	R = corrcoef(yLag(npMask), yNext(npMask));
	out.ac1NP = R(2, 1);
else
	out.ac1NP = NaN;
end
out.ac1AsymCross = out.ac1PN - out.ac1NP;

end
