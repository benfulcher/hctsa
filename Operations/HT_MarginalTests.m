function p = HT_MarginalTests(y, theTest)
% HT_MarginalTests     p-value of a hypothesis test about the marginal distribution of values.
%
% These tests ask a question about the distribution of the values in the time
% series (its center, spread, symmetry, or shape) and are insensitive to the
% temporal ordering of the measurements: any reordering of the series gives the
% same p-value. For tests of the temporal (serial) structure of the series, see
% HT_IndependenceTests.
%
% Tests are implemented as functions in Matlab's Statistics Toolbox.
%
% ---INPUTS:
% y, the input time series
%
% theTest, the hypothesis test to perform (and the null hypothesis that it tests):
%           (i) 'signtest': sign test, the data are a random sample from a
%                       continuous distribution with a median of zero
%           (ii) 'signrank': Wilcoxon signed rank test, the data are a random sample
%                       from a continuous, symmetric distribution with a median of zero
%           (iii) 'jbtest': Jarque-Bera test of composite normality, the data are
%                       drawn from a normal distribution with unknown mean and variance
%                       (the test statistic is based on the sample skewness and kurtosis)
%           (iv) 'vartest': variance test, the data are drawn from a normal
%                       distribution with a variance of one (and unknown mean)
%           (v) 'ztest': Z-test, the data are drawn from a normal distribution with
%                       a mean of zero and a (known) standard deviation of one
%
% ---OUTPUT:
% The p-value of the specified test: the probability, under the null hypothesis
% above, of a test statistic at least as extreme as that observed. Small values
% are evidence against the null hypothesis.

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
	case 'signtest' % Statistics Toolbox
		[p, ~] = signtest(y);
		% for some reason this one has p-value as the first output

	case 'vartest' % Statistics Toolbox
		[~, p] = vartest(y, 1); % normal distribution of variance 1

	case 'ztest' % Statistics Toolbox
		[~, p] = ztest(y, 0, 1);

	case 'signrank' % Statistics Toolbox
		[p, ~] = signrank(y);

	case 'jbtest' % Statistics Toolbox
		warning('off', 'stats:jbtest:PTooBig'); % suspend this warning
		warning('off', 'stats:jbtest:PTooSmall'); % suspend this warning
		[~, p] = jbtest(y);
		warning('on', 'stats:jbtest:PTooBig'); % resume this warning
		warning('on', 'stats:jbtest:PTooSmall'); % resume this warning

	otherwise
		error('Unknown hypothesis test ''%s''', theTest);
end

end
