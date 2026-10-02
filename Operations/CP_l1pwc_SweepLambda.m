function out = CP_l1pwc_SweepLambda(y, lambdar)
% CP_l1pwc_SweepLambda   How step detection changes as its penalty on step size varies.
%
% Gives information about discrete steps in the signal across a range of
% regularization parameters lambda, using the l1pwc function from Max Little's
% step-detection toolkit. At each lambda in lambdar, CP_ML_StepDetect(y, 'l1pwc',
% lambda) is run, and the number of constant segments per sample (nsegments), the
% reduction in standard deviation from removing the piecewise-constant fit
% (rmsoff), and that reduction per constant segment (rmsoffpstep) are recorded.
% The outputs summarize how these quantities vary with lambda. Note that a
% lambda below 1 is taken by CP_ML_StepDetect as a proportion of the largest
% lambda that gives any steps.
%
% ---INPUTS:
% y, the input time series
%
% lambdar, a vector specifying the lambda parameters to use
%
% ---OUTPUTS:
% rmserrsu05, rmserrsu02, rmserrsu01: the first lambda in lambdar at which the
%       reduction in standard deviation (rmsoff) falls below 0.5, 0.2, 0.1
% nsegsu005, nsegsu001: the first lambda in lambdar at which the number of
%       segments per sample falls below 0.05, 0.01
%       (all five are NaN if the threshold is never crossed)
% corrsegerr, the correlation across lambdar between the number of segments and
%       the reduction in standard deviation
% bestrmserrpseg, the maximum reduction in standard deviation per constant
%       segment (rmsoffpstep) over lambdar
% bestlambda, the lambda at which that maximum occurs
%
% ---REFERENCES:
% Max A. Little and Nick S. Jones, "Sparse Bayesian Step-Filtering for
% High-Throughput Analysis of Molecular Machine Dynamics", Proc. ICASSP (2010).

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

Llambdar = length(lambdar);
nsegs = zeros(Llambdar, 1);
rmserrs = zeros(Llambdar, 1);
rmserrpsegs = zeros(Llambdar, 1);

for i = 1:length(lambdar)
	lambda = lambdar(i);
	% Run the (stochastic) step detection algorithm:
	outi = CP_ML_StepDetect(y, 'l1pwc', lambda);
	nsegs(i) = outi.nsegments;
	rmserrs(i) = outi.rmsoff;
	rmserrpsegs(i) = outi.rmsoffpstep;
end

% ------------------------------------------------------------------------------
% Define the FirstUnder function:
% ------------------------------------------------------------------------------
% Often finding the first time a certain vector (x) drops under a given
% threshold (y). This function does it:

FirstUnder = @(x, y) find(x < y, 1, 'first');

% ------------------------------------------------------------------------------
% (*) Use rmsunderx subfunction to analyze when RMS errors drop under a set of
% thresholds, *x*, for the first time:
% ------------------------------------------------------------------------------
out.rmserrsu05 = NaNIfEmpty(lambdar(FirstUnder(rmserrs, 0.5)));
out.rmserrsu02 = NaNIfEmpty(lambdar(FirstUnder(rmserrs, 0.2)));
out.rmserrsu01 = NaNIfEmpty(lambdar(FirstUnder(rmserrs, 0.1)));

% ------------------------------------------------------------------------------
% (*) Use the nsegsunderx subfunction to analyze when nseg drops under a set of
% thresholds, *x*, for the first time:
% nsegunderx = @(x) find(nsegs < x, 1, 'first');
% ------------------------------------------------------------------------------
out.nsegsu005 = NaNIfEmpty(lambdar(FirstUnder(nsegs, 0.05)));
out.nsegsu001 = NaNIfEmpty(lambdar(FirstUnder(nsegs, 0.01)));

% Calculate the correlation between the number of segments and rmserrs
R = corrcoef(nsegs, rmserrs);
out.corrsegerr = R(2, 1);

% Maximum rmserrpsegment
indbest = find(rmserrpsegs == max(rmserrpsegs), 1, 'first'); % where the maximum occurs
out.bestrmserrpseg = rmserrpsegs(indbest);
out.bestlambda = lambdar(indbest);

% ------------------------------------------------------------------------------
% Subfunctions
% ------------------------------------------------------------------------------
function N = NaNIfEmpty(x)
	% Returns a NaN if x is empty, otherwise returns x.
	if isempty(x)
		N = NaN;
	else
		N = x;
	end
end

end
