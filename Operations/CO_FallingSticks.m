function out = CO_FallingSticks(y)
% CO_FallingSticks   Physical falling-sticks model of line-of-sight interaction.
%
% As in CO_StickAngles, each time-series value is treated as a rigid stick
% standing on the zero baseline, with sticks grouped by sign into a
% 'positive' set (protruding up from the zero level) and a 'negative' set
% (protruding down). Here, sticks are toppled: each stick rotates about its
% base towards later same-sign sticks and stops at whichever angle first
% brings it into contact with one -- either its trunk striking the side of
% a taller later stick, or its underside striking the tip of a shorter one
% it topples clean over -- or else it falls flat (angle = pi/2) if no
% later same-sign stick lies within reach (a stick of height h can only
% ever reach as far as horizontal distance h).
%
% This differs from CO_StickAngles, which only ever compares a stick to
% its immediate same-sign successor via the slope between them.
% CO_FallingSticks instead allows a stick to skip over intervening sticks
% to hit a farther one, so it is sensitive to range-dependent local-
% extremum structure (e.g. a tall stick toppling clean over an
% intervening short one) that the purely local (i, i+1) comparison cannot
% see.
%
% Adapted from a Python 'FALLstick' reference implementation by Eugene Chon
% <eugenechon04@gmail.com>.
%
% ---INPUTS:
% y, the input time series (assumed z-scored: the sign split is around the
%       mean, matching CO_StickAngles's convention; heights are in z-scored units
%       and horizontal distances in samples)
%
% ---OUTPUTS:
% Fall angles are in radians from the vertical: 0 is upright, pi/2 is flat. The last
% stick of each sign has no later stick and contributes no angle (but counts as
% falling flat). Suffixes _p, _n and _all refer to the positive sticks, the
% negative sticks and both.
% mean_p, median_p, mean_n, median_n, mean_all, median_all, the mean and median fall
%       angle,
% std_all, the standard deviation of all fall angles,
% diff_pn, mean_p - mean_n, the asymmetry between positive and negative sticks,
% propFlat_p, propFlat_n, propFlat_all, the proportion of sticks that fall flat,
% propSkip_p, propSkip_n, propSkip_all, of the sticks that hit another stick, the
%       proportion that skip over one or more sticks to hit a farther one,
% propCase2_p, propCase2_n, propCase2_all, the proportion of sticks whose fall ends
%       in a tip-strike (topple-over, case 2) rather than a trunk-strike (case 1);
%       flat falls count as case 1,
% skewness_all, kurtosis_all, q10_all, the skewness, kurtosis and 10th percentile
%       of all fall angles (the 90th percentile is omitted: it was always the
%       pi/2 flat-fall spike),
% tau_p, tau_n, the first zero-crossing of the autocorrelation function of the
%       sequence of fall angles of the positive or negative sticks (interpolated),
% ac1_p, ac1_n, the lag-1 autocorrelation of those sequences (tau and ac1 are NaN
%       when the angles are constant, to rounding error, e.g. when all sticks fall flat;
%       likewise skewness_all and kurtosis_all, and std_all is then 0).
% Statistics with no sticks to work on are NaN.

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

ixPos = find(y >= 0);
ixNeg = find(y < 0);

[anglesPos, colourPos, casePos] = SUB_fallBranch(ixPos, y);
[anglesNeg, colourNeg, caseNeg] = SUB_fallBranch(ixNeg, y);

allAngles = [anglesPos; anglesNeg];

% ------------------------------------------------------------------------------
%% Location and spread of the fall-angle distribution
% ------------------------------------------------------------------------------
out.mean_p = SUB_safeStat(@mean, anglesPos);
out.median_p = SUB_safeStat(@median, anglesPos);
out.mean_n = SUB_safeStat(@mean, anglesNeg);
out.median_n = SUB_safeStat(@median, anglesNeg);

out.mean_all = SUB_safeStat(@mean, allAngles);
out.median_all = SUB_safeStat(@median, allAngles);
out.std_all = SUB_safeStat(@std, allAngles);
if ~isempty(allAngles) && SUB_isConstant(allAngles)
	out.std_all = 0; % constant to rounding error: the spread is zero
end

% Asymmetry between the positive- and negative-branch fall angles:
if ~isnan(out.mean_p) && ~isnan(out.mean_n)
	out.diff_pn = out.mean_p - out.mean_n;
else
	out.diff_pn = NaN;
end

% ------------------------------------------------------------------------------
%% Collision-type proportions: flat / immediate-hit / skip-hit
% ------------------------------------------------------------------------------
out.propFlat_p = SUB_propFirst(colourPos);
out.propFlat_n = SUB_propFirst(colourNeg);
out.propFlat_all = SUB_propFirstPooled(colourPos, colourNeg);

out.propSkip_p = SUB_propSkip(colourPos);
out.propSkip_n = SUB_propSkip(colourNeg);
out.propSkip_all = SUB_propSkipPooled(colourPos, colourNeg);

% ------------------------------------------------------------------------------
%% Hit-type proportions: trunk-strike (case 1) vs. tip-strike/topple-over (case 2)
% ------------------------------------------------------------------------------
out.propCase2_p = SUB_propCase2(casePos);
out.propCase2_n = SUB_propCase2(caseNeg);
out.propCase2_all = SUB_propCase2Pooled(casePos, caseNeg);

% ------------------------------------------------------------------------------
%% Distribution shape
% ------------------------------------------------------------------------------
% (The 90th percentile is omitted: across two independent redundancy-check
% datasets it landed exactly on the pi/2 flat-fall spike every time, since
% most series have >=10% flat falls -- it carries no information.)
if length(allAngles) >= 2
	out.q10_all = quantile(allAngles, 0.1);
else
	out.q10_all = NaN;
end
if length(allAngles) >= 2 && ~SUB_isConstant(allAngles)
	out.skewness_all = skewness(allAngles);
	out.kurtosis_all = kurtosis(allAngles);
else
	out.skewness_all = NaN; % undefined for constant angles
	out.kurtosis_all = NaN;
end

% ------------------------------------------------------------------------------
%% Persistence of the fall-angle sequence
% ------------------------------------------------------------------------------
if length(anglesPos) >= 2 && ~SUB_isConstant(anglesPos)
	zAnglesPos = zscore(anglesPos);
	out.tau_p = CO_FirstCrossing(zAnglesPos, 'ac', 0, 'continuous');
	out.ac1_p = CO_AutoCorr(zAnglesPos, 1, 'Fourier');
else
	out.tau_p = NaN;
	out.ac1_p = NaN;
end

if length(anglesNeg) >= 2 && ~SUB_isConstant(anglesNeg)
	zAnglesNeg = zscore(anglesNeg);
	out.tau_n = CO_FirstCrossing(zAnglesNeg, 'ac', 0, 'continuous');
	out.ac1_n = CO_AutoCorr(zAnglesNeg, 1, 'Fourier');
else
	out.tau_n = NaN;
	out.ac1_n = NaN;
end

% -------------------------------------------------------------------------------
function isConst = SUB_isConstant(x)
	% A sequence of angles counts as constant when its spread is at rounding level
	% (e.g., all sticks fall flat, but some angles differ in the last bits): zscore
	% would then amplify the rounding noise, so the persistence statistics are NaN.
	isConst = std(x) <= 1e-10 * max(1, max(abs(x)));
end

function [angles, colourCounts, caseCounts] = SUB_fallBranch(ix, y)
	% Topples each stick in a same-sign index list ix towards later sticks in
	% the same list, and returns:
	%   angles       -- fall angle for each processed stick (all but the last,
	%                    which trivially always falls flat and so contributes
	%                    no angle of its own)
	%   colourCounts -- [numFlat, numImmediateHit, numSkipHit], summing to
	%                    length(ix) (the forced-flat last stick included)
	%   caseCounts   -- [numCase1, numCase2], summing to length(ix) (the
	%                    forced-flat last stick counted as case 1)
	nj = length(ix);
	if nj == 0
		% No sticks at all in this branch (e.g. an all-positive or
		% all-negative series): nothing falls, so there is no phantom
		% 'last bar' to count as a flat fall either.
		angles = zeros(0, 1);
		colourCounts = [0, 0, 0];
		caseCounts = [0, 0];
		return
	end
	if nj == 1
		angles = zeros(0, 1);
		colourCounts = [1, 0, 0]; % the lone stick always falls flat
		caseCounts = [1, 0];
		return
	end

	angles = zeros(nj - 1, 1);
	colourFlag = zeros(nj - 1, 1); % 0 = flat, 1 = immediate hit, 2 = skip hit
	caseFlag = zeros(nj - 1, 1); % 1 or 2

	for i = 1:nj - 1
		x1 = ix(i);
		y1 = y(x1);
		height1 = abs(y1);

		minAngle = pi/2; % default: falls flat
		fallK = i; % lands on itself if flat
		minCase1 = true;

		if height1 > 0 % a zero-height stick is already lying flat
			for k = i + 1:nj
				x2 = ix(k);
				y2 = y(x2);
				dx = x2 - x1; % a stick of height1 can reach at most height1
				if dx > height1
					break % all later sticks are farther still -- out of reach
				end

				% Case 1 (trunk-strike): stick 1's tip hits the side of stick
				% 2. Case 2 (tip-strike): stick 1 is taller than stick 2 and
				% would clear its top before reaching dx, so instead its
				% underside comes down onto stick 2's tip.
				isCase1 = true;
				if abs(y1) > abs(y2) && height1 * sin(acos(y2 / y1)) > dx
					isCase1 = false;
				end

				if isCase1
					angle = asin(dx / height1);
				else
					angle = atan(dx / abs(y2));
				end
				angle = min(angle, pi - angle);

				if angle < minAngle
					minAngle = angle;
					fallK = k;
					minCase1 = isCase1;
				end
			end
		end

		angles(i) = minAngle;
		if fallK == i
			colourFlag(i) = 0;
		elseif fallK == i + 1
			colourFlag(i) = 1;
		else
			colourFlag(i) = 2;
		end
		caseFlag(i) = 2 - minCase1; % true (case 1) -> 1, false (case 2) -> 2
	end

	colourCounts = [sum(colourFlag == 0) + 1, sum(colourFlag == 1), sum(colourFlag == 2)];
	caseCounts = [sum(caseFlag == 1) + 1, sum(caseFlag == 2)];
end

% -------------------------------------------------------------------------------
function statOut = SUB_safeStat(f, x)
	if isempty(x)
		statOut = NaN;
	else
		statOut = f(x);
	end
end

function p = SUB_propFirst(colourCounts)
	total = sum(colourCounts);
	if total == 0
		p = NaN;
	else
		p = colourCounts(1) / total;
	end
end

function p = SUB_propFirstPooled(colourA, colourB)
	total = sum(colourA) + sum(colourB);
	if total == 0
		p = NaN;
	else
		p = (colourA(1) + colourB(1)) / total;
	end
end

function p = SUB_propSkip(colourCounts)
	nonFlat = colourCounts(2) + colourCounts(3);
	if nonFlat == 0
		p = NaN;
	else
		p = colourCounts(3) / nonFlat;
	end
end

function p = SUB_propSkipPooled(colourA, colourB)
	nonFlat = colourA(2) + colourA(3) + colourB(2) + colourB(3);
	if nonFlat == 0
		p = NaN;
	else
		p = (colourA(3) + colourB(3)) / nonFlat;
	end
end

function p = SUB_propCase2(caseCounts)
	total = sum(caseCounts);
	if total == 0
		p = NaN;
	else
		p = caseCounts(2) / total;
	end
end

function p = SUB_propCase2Pooled(caseA, caseB)
	total = sum(caseA) + sum(caseB);
	if total == 0
		p = NaN;
	else
		p = (caseA(2) + caseB(2)) / total;
	end
end

end
