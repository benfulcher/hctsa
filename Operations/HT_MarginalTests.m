function p = HT_MarginalTests(y, theTest)
% HT_MarginalTests   p-value of a hypothesis test about the distribution of values.
%
% These tests ask a question about the distribution of the values in the time
% series (its center, spread, symmetry, or shape) and are insensitive to the
% temporal ordering of the measurements: any reordering of the series gives the
% same p-value. For tests of the temporal (serial) structure of the series, see
% HT_IndependenceTests.
%
% Tests are implemented as functions in MATLAB's Statistics Toolbox.
%
% ---INPUTS:
% y, the input time series
% theTest, the hypothesis test to perform (and the null hypothesis that it
%       tests):
%           (i) 'signtest': sign test, the data are a random sample from a
%                       continuous distribution with a median of zero
%           (ii) 'signrank': Wilcoxon signed rank test, the data are a random
%                       sample from a continuous, symmetric distribution with
%                       a median of zero
%           (iii) 'jbtest': Jarque-Bera test of composite normality, the data
%                       are drawn from a normal distribution with unknown mean
%                       and variance (the test statistic is based on the
%                       sample skewness and kurtosis)
%           (iv) 'vartest': variance test, the data are drawn from a normal
%                       distribution with a variance of one (and unknown mean)
%           (v) 'ztest': Z-test, the data are drawn from a normal distribution
%                       with a mean of zero and a (known) standard deviation of
%                       one
%           (vi) 'logjb': the logarithm of the Jarque-Bera statistic,
%                       JB = n/6*(S^2 + (K - 3)^2/4), with S and K the (biased)
%                       sample skewness and kurtosis (Jarque and Bera 1987): how far
%                       the shape of the distribution is from normal. Unlike the
%                       'jbtest' p-value, which is read from a table of limited range
%                       (clamped to [0.001, 0.5]), it is not bounded, so it
%                       distinguishes strongly non-normal distributions. NaN for a
%                       constant series (S and K undefined) or if JB is exactly 0.
%
% ---OUTPUTS:
% p, the p-value of the specified test: the probability, under the null
% hypothesis above, of a test statistic at least as extreme as that observed
% (for 'logjb', log(JB) instead).
%
% ---REFERENCES:
% C. M. Jarque and A. K. Bera, "A test for normality of observations and
% regression residuals", International Statistical Review 55(2), 163-172 (1987).
% DOI: 10.2307/1403192
% Small values are evidence against the null hypothesis.

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

	case 'logjb'
		y = y(:);
		n = length(y);
		d = y - mean(y);
		m2 = mean(d.^2);
		S = mean(d.^3) / m2^1.5; % sample skewness
		K = mean(d.^4) / m2^2; % sample kurtosis
		JB = n / 6 * (S^2 + (K - 3)^2 / 4);
		if JB > 0 && isfinite(JB)
			p = log(JB);
		else
			p = NaN; % constant series, or exactly normal skewness and kurtosis (log(0))
		end

	otherwise
		error('Unknown hypothesis test ''%s''', theTest);
end

end
