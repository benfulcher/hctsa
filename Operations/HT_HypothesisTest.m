function p = HT_HypothesisTest(x, theTest)
% HT_HypothesisTest     [DEPRECATED] Statistical hypothesis test applied to a time series.
%
% DEPRECATED: use HT_MarginalTests (tests about the distribution of values)
% or HT_IndependenceTests (tests of serial independence) instead; this thin
% wrapper is kept only so that custom input files keep working, and is no longer
% part of the default feature library.
%
% ---INPUTS:
% x, the input time series
%
% theTest, the hypothesis test to perform, dispatched as:
%           HT_MarginalTests: 'signtest', 'vartest', 'ztest', 'signrank', 'jbtest'
%           HT_IndependenceTests: 'runstest', 'lbq'
%
% ---OUTPUT:
% p-value from the specified statistical test (identical to that of the function
% it dispatches to)

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

% Warn once per session:
persistent hasWarned
if isempty(hasWarned)
	warning('hctsa:deprecated', ['HT_HypothesisTest is deprecated: use HT_MarginalTests ' ...
				'(distribution tests) or HT_IndependenceTests (serial independence tests) instead.']);
	hasWarned = true;
end

switch theTest
	case {'signtest','vartest','ztest','signrank','jbtest'}
		p = HT_MarginalTests(x, theTest);

	case {'runstest','lbq'}
		p = HT_IndependenceTests(x, theTest);

	otherwise
		error('Unknown hypothesis test ''%s''', theTest);
end

end
