function out = TSTL_delaytime(y, maxDelay, past, randomSeed)
% TSTL_delaytime    Optimal delay time using the method of Parlitz and Wichard.
%
% Computed natively in MATLAB. TSTOOL's own 'delaytime' code
% (tstoolbox/@signal/delaytime.m, vendored in this repo under
% Toolboxes/OpenTSTOOL) turned out to be pure MATLAB itself -- no compiled
% TSTOOL binary/mex involved at all -- so its exact algorithm (still
% undocumented/uncredited beyond "method of Parlitz and Wichard", per
% TSTOOL's own comments) is reproduced here directly, rather than
% approximated: sort the (truncated) series by value; repeatedly (64
% iterations) pick a random reference point by value-rank, then find its
% nearest value-neighbors on either side (excluding a Theiler window of
% "past" samples in time) and accumulate |x(neighbor+lag) - x(ref+lag)|
% for lag = 0:maxDelay; the average over iterations is tau(lag+1).
%
% Two bugs fixed relative to TSTOOL's own source's search for the
% nearest-by-value neighbor *above* the reference rank
% ("post = index(find(abs(index(ref+1:end)-actual) > past))"), confirmed
% directly with a concrete worked example before fixing:
% (1) it applies find() to the subarray index(ref+1:end) but then indexes
%     back into the *full* index array with the resulting (subarray-
%     relative) positions without adding ref back on -- so it actually
%     resamples low ranks (1:length(index(ref+1:end))), the same region
%     "pre" already searches, rather than the intended ranks above ref;
% (2) even with that offset fixed, taking the last element of the
%     (ascending-rank) candidates above ref gives the *farthest* passing
%     neighbor, not the nearest -- the mirror image of how "pre" already
%     takes its last (ascending-rank, i.e. nearest-to-ref) candidate below
%     ref. "post" here takes the first element after offsetting, so both
%     directions consistently return the nearest value-neighbor passing
%     the Theiler-window exclusion.
%
% ---INPUTS:
% y, column vector of time series data
%
% maxDelay, maximum value of the delay to consider (can also specify a
%           proportion of time series length)
%
% past, the TSTOOL documentation describes this parameter as "?", which is
%       relatively uninformative.
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
%% Preliminaries
% ------------------------------------------------------------------------------
y = y(:); % ensure a column vector
N = length(y); % length of time series

% ------------------------------------------------------------------------------
% Check Inputs:
% ------------------------------------------------------------------------------
% (1) Maximum delay, maxDelay
if nargin < 2 || isempty(maxDelay)
	maxDelay = 0.2; % 1/5 the length of the time series
end
if maxDelay < 1 && maxDelay > 0
	maxDelay = round(N * maxDelay); % specify a proportion of time series length
end

if maxDelay < 10
	maxDelay = 10;
	fprintf(1, 'Max delay set to its minimum: delaytime = 10\n');
end
if maxDelay >= N / 2
	% Heuristic for appropriate time delay
	warning('Max delay, %u, too long for time series of length %u', maxDelay, N)
	out = NaN; return
end

% randomSeed: how to treat the randomization
if nargin < 4
	randomSeed = [];
end

% ------------------------------------------------------------------------------
%% Run
% ------------------------------------------------------------------------------

% Control the random seed (for reproducibility):
BF_ResetSeed(randomSeed);

% Reproduces TSTOOL's tstoolbox/@signal/delaytime.m directly (see header comment):
ITERATIONS = 64;
len = N - maxDelay;
[~, index] = sort(y(1:len));

err = zeros(maxDelay + 1, 1);
for i = 1:ITERATIONS
	while true
		ref = ceil(rand(1, 1) * len);
		actual = index(ref);
		preCandidates = index(abs(index(1:ref - 1) - actual) > past);
		postCandidates = index(ref + find(abs(index(ref + 1:end) - actual) > past));
		if ~isempty(preCandidates) && ~isempty(postCandidates)
			pre = preCandidates(end); % nearest-in-value candidate below ref
			post = postCandidates(1); % nearest-in-value candidate above ref
			break
		end
	end
	errpre = abs(y(pre:pre + maxDelay) - y(actual:actual + maxDelay));
	errpost = abs(y(post:post + maxDelay) - y(actual:actual + maxDelay));
	err = err + errpre + errpost;
end

tau = err / ITERATIONS;

% plot(tau);

% ------------------------------------------------------------------------------
%% Output Statistics
% ------------------------------------------------------------------------------
% tau tends to start low and then rise to some (noisy) value
out.tau1 = tau(1);
out.tau2 = tau(2);
out.tau3 = tau(3);
out.difftau12 = tau(2) - tau(1);
out.difftau13 = tau(3) - tau(1);
out.meantau = mean(tau);
out.stdtau = std(tau);
out.mintau = min(tau);
out.maxtau = max(tau);

end
