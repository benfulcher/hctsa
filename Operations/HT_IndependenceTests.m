function p = HT_IndependenceTests(y, theTest)
% HT_IndependenceTests     p-value of a hypothesis test of serial independence.
%
% These tests ask whether the temporal ordering of the time series carries
% structure: the null hypothesis is that successive values are independent (no
% serial dependence). Unlike the tests in HT_MarginalTests, the p-value
% depends on the order of the measurements, and can change if the series is
% reordered.
%
% Tests are implemented as functions in Matlab's Statistics Toolbox
% (except the Ljung-Box Q-test, which uses the Econometrics Toolbox).
%
% ---INPUTS:
% y, the input time series
%
% theTest, the hypothesis test to perform (and the null hypothesis that it tests):
%           (i) 'runstest': runs test for randomness, the values occur in random
%                       order (assessed from the number of runs of values above and
%                       below a cutoff, using the Matlab default cutoff, the median)
%           (ii) 'lbq': Ljung-Box Q-test for residual autocorrelation, the series
%                       has no autocorrelation (all autocorrelations are zero,
%                       jointly over the lags considered by Matlab's default)
%
% ---OUTPUTS:
% p, a scalar: the p-value of the specified test: the probability, under the null hypothesis
% above, of a test statistic at least as extreme as that observed. Small values
% are evidence against the null hypothesis (i.e., evidence of serial dependence).

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

switch theTest
	case 'runstest' % Statistics Toolbox
		[~, p] = runstest(y);

	case 'lbq'
		% Check that an Econometrics Toolbox license is available:
		BF_CheckToolbox('econometrics_toolbox');

		% Perform the test
		[~, p] = lbqtest(y);

	otherwise
		error('Unknown hypothesis test ''%s''', theTest);
end

end
