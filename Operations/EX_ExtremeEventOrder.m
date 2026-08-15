function out = EX_ExtremeEventOrder(y, extremeThresh)
% EX_ExtremeEventOrder   Temporal patterning of positive vs. negative extreme events.
%
% Labels each point in the top/bottom extremeThresh fraction of the
% distribution as a positive- or negative-direction 'extreme event', then
% asks not just how often each type occurs (that's DN_OutlierInclude's
% territory), but how the two types are *ordered* in time relative to each
% other: does a positive extreme tend to be followed by another positive
% one (clustering by sign), or does the sequence alternate (P,N,N,P,N,...)?
%
% ---INPUTS:
% y, the input time series (assumed z-scored)
%
% extremeThresh, the proportion of points (in each direction) to count as
%       'extreme' (default: 0.05, i.e., the most extreme 5% each way)
%
% ---OUTPUTS:
% propPosEvents: what proportion of extreme events are positive-direction.
%       (the total event count itself is not reported here: with events
%       defined by a fixed quantile threshold, numEvents is just
%       extremeThresh*2*N up to tie noise -- a near-exact re-encoding of
%       series length, not a property of the dynamics.)
% meanInterval, cvInterval: mean and coefficient of variation of the
%       inter-event intervals of the combined (both-direction) event
%       sequence -- how bursty extreme events are overall, irrespective of
%       direction.
% meanIntervalPos, cvIntervalPos, meanIntervalNeg, cvIntervalNeg: the same,
%       but computed separately within each direction's own event
%       sub-sequence (e.g. meanIntervalPos is the mean gap between
%       consecutive positive-direction events, ignoring any negative
%       events that occur in between) -- this is a different question from
%       the combined interval stats above: a series can have very bursty
%       positive extremes and very regular negative extremes even while
%       the combined (direction-blind) event stream looks unremarkable.
% alternationRate: proportion of consecutive event pairs whose direction
%       differs -- 1 means the sequence strictly alternates sign, 0 means
%       same-direction events always cluster together.
% propPN, propNP: alternationRate split by switch direction -- the
%       proportion of consecutive event pairs that switch positive-to-
%       negative, and negative-to-positive, respectively (propPN+propNP ==
%       alternationRate). Mirrors CO_PosNegAsymmetry's propPN/propNP.
% meanIntervalPN, cvIntervalPN, meanIntervalNP, cvIntervalNP: mean and
%       coefficient of variation of the gaps between successive PN
%       switches (and, separately, successive NP switches), timestamped at
%       the later (post-switch) event. Different question again from
%       meanIntervalPos/Neg above: those count every event of a given
%       direction; these count only the (rarer) moments the direction
%       actually flips, e.g. a series alternating in bursts P,P,P,N,N,N,...
%       has short meanIntervalPos/Neg (events of a given type are close
%       together within a burst) but long meanIntervalPN/NP (switches
%       between bursts are rare).
% lzComplexity: normalized Lempel-Ziv complexity (via Michael Small's
%       MS_complexitybs, called directly on the true 0/1 label sequence --
%       NOT via MS_complexity/EN_LZComplexity, whose equiprobable
%       re-binning would relabel points to force equal counts per symbol,
%       corrupting the true P/N assignment whenever propPosEvents differs
%       from 0.5) of the event-direction sequence: captures patterning
%       beyond pairwise alternation alone, e.g. a period-4 repeat
%       P,P,N,N,P,P,N,N,... has the same alternationRate (~0.5) as i.i.d.
%       P/N noise, but is far more predictable (much lower lzComplexity).

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

if nargin < 2 || isempty(extremeThresh)
	extremeThresh = 0.05; % most extreme 5% of points, each direction
end
if extremeThresh <= 0 || extremeThresh >= 0.5
	error('extremeThresh must be in (0,0.5) so the two tails cannot overlap.');
end

y = y(:);

% ------------------------------------------------------------------------------
% Identify extreme events by direction (raw threshold exceedances, not
% declustered peaks -- runs of consecutive same-direction exceedances are
% themselves part of the signal this operation is trying to characterize,
% e.g. persistent positive-extreme regimes should show up as low
% alternationRate, not be collapsed away by declustering).
% ------------------------------------------------------------------------------
upperThresh = quantile(y, 1 - extremeThresh);
lowerThresh = quantile(y, extremeThresh);

posIdx = find(y > upperThresh);
negIdx = find(y < lowerThresh);

eventIdx = [posIdx; negIdx];
eventType = [ones(length(posIdx), 1); zeros(length(negIdx), 1)]; % 1 = positive, 0 = negative
[eventIdx, sortOrder] = sort(eventIdx);
eventType = eventType(sortOrder);

numEvents = length(eventIdx);
out.propPosEvents = mean(eventType);

% ------------------------------------------------------------------------------
% Event timing: combined (both directions), and separately within each
% direction's own sub-sequence of events (posIdx/negIdx are already sorted,
% since find() returns indices in increasing order).
% ------------------------------------------------------------------------------
[out.meanInterval, out.cvInterval] = SUB_intervalStats(eventIdx);
[out.meanIntervalPos, out.cvIntervalPos] = SUB_intervalStats(posIdx);
[out.meanIntervalNeg, out.cvIntervalNeg] = SUB_intervalStats(negIdx);

% ------------------------------------------------------------------------------
% Ordering of event direction: alternation rate and Lempel-Ziv complexity
% ------------------------------------------------------------------------------
if numEvents >= 2
	typeChange = diff(eventType); % +1 = N-to-P switch, -1 = P-to-N switch, 0 = no switch
	out.alternationRate = mean(typeChange ~= 0);
	out.propPN = mean(typeChange == -1);
	out.propNP = mean(typeChange == 1);
	% Timestamp each switch at the later (post-switch) event:
	pnSwitchIdx = eventIdx(find(typeChange == -1) + 1);
	npSwitchIdx = eventIdx(find(typeChange == 1) + 1);
else
	out.alternationRate = NaN;
	out.propPN = NaN;
	out.propNP = NaN;
	pnSwitchIdx = [];
	npSwitchIdx = [];
end
[out.meanIntervalPN, out.cvIntervalPN] = SUB_intervalStats(pnSwitchIdx);
[out.meanIntervalNP, out.cvIntervalNP] = SUB_intervalStats(npSwitchIdx);

% Below ~10 symbols, normalized LZ complexity (itself an asymptotic
% comparison against the expected count for a noise sequence) becomes
% unstable, so leave it NaN rather than report a number with no meaningful
% floor to compare against. Note this is a bare minimum, not a guarantee of
% an unbiased estimate -- like other symbolic-complexity measures, it
% retains a mild downward bias for short sequences that only settles by a
% few hundred events; interpret cautiously for short input series.
%
% Called with a single argument: MS_complexitybs's compiled mex silently
% ignores its documented second ('n', alphabet size) argument and instead
% infers the alphabet size from the data itself (floor(x)+1, then the max
% over the sequence) -- correct here since eventType is genuinely 0/1, but
% passing 'n' just triggers a spurious "Second input ignored" warning.
if numEvents >= 10 && out.propPosEvents > 0 && out.propPosEvents < 1
	out.lzComplexity = MS_complexitybs(eventType);
else
	% propPosEvents in {0,1} means every event has the same direction (a
	% single symbol value) -- MS_complexitybs's normalization divides by
	% log(numDistinctSymbols), which is log(1)=0 in that case, so it would
	% return Inf rather than a meaningful complexity value.
	out.lzComplexity = NaN;
end

end

% ------------------------------------------------------------------------------
function [meanI, cvI] = SUB_intervalStats(idx)
	% Mean and coefficient of variation of the gaps between consecutive
	% (already time-sorted) event indices in idx.
	if length(idx) >= 2
		intervals = diff(idx);
		meanI = mean(intervals);
		if meanI > 0
			cvI = std(intervals) / meanI;
		else
			cvI = NaN;
		end
	else
		meanI = NaN;
		cvI = NaN;
	end
end
