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
% numEvents, propPosEvents: how many extreme events there are in total, and
%       what proportion of them are positive-direction.
% meanInterval, cvInterval: mean and coefficient of variation of the
%       inter-event intervals of the combined (both-direction) event
%       sequence -- how bursty extreme events are overall, irrespective of
%       direction.
% alternationRate: proportion of consecutive event pairs whose direction
%       differs -- 1 means the sequence strictly alternates sign, 0 means
%       same-direction events always cluster together.
% lzComplexity: normalized Lempel-Ziv complexity (via Michael Small's
%       MS_complexitybs, called directly on the true 0/1 label sequence --
%       NOT via MS_complexity/EN_MS_LZcomplexity, whose equiprobable
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
out.numEvents = numEvents;
out.propPosEvents = mean(eventType);

% ------------------------------------------------------------------------------
% Overall event timing (both directions combined)
% ------------------------------------------------------------------------------
if numEvents >= 2
	intervals = diff(eventIdx);
	out.meanInterval = mean(intervals);
	if out.meanInterval > 0
		out.cvInterval = std(intervals) / out.meanInterval;
	else
		out.cvInterval = NaN;
	end
else
	out.meanInterval = NaN;
	out.cvInterval = NaN;
end

% ------------------------------------------------------------------------------
% Ordering of event direction: alternation rate and Lempel-Ziv complexity
% ------------------------------------------------------------------------------
if numEvents >= 2
	out.alternationRate = mean(diff(eventType) ~= 0);
else
	out.alternationRate = NaN;
end

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
