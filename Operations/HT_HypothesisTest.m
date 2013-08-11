% HT_HypothesisTest
% 
% Outputs the p-value from a statistical hypothesis test applied to the
% time series.
% 
% Tests are implemented as functions in Matlab's Statistics Toolbox.
% (except Ljung-Box Q-test, which uses the Econometrics Toolbox)
% 
% INPUTS:
% x, the input time series
% 
% thetest, the hypothesis test to perform:
%           (i) sign test ('signtest'),
%           (ii) runs test ('runstest'),
%           (iii) variance test ('vartest'),
%           (iv) Z-test ('ztest'),
%           (v) Wilcoxon signed rank test for a zero median ('signrank'),
%           (vi) Jarque-Bera test of composite normality ('jbtest').
%           (vii) Ljung-Box Q-test for residual autocorrelation ('lbq')
%           
% OUTPUT:
% the p-value from the statistical test
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function p = HT_HypothesisTest(x,thetest)
% Ben Fulcher, 2009

switch thetest
    case 'signtest' % Statistics Toolbox
        [p, ~] = signtest(x);
        % for some reason this one has p-value as the first output
        
    case 'runstest' % Statistics Toolbox
        [~, p] = runstest(x);

    case 'vartest' % Statistics Toolbox
        [~, p] = vartest(x,1); % normal distribution of variance 1
        
    case 'ztest' % Statistics Toolbox
        [~, p] = ztest(x,0,1);
        
    case 'signrank' % Statistics Toolbox
        [p, ~] = signrank(x);
        
    case 'jbtest' % Statistics Toolbox
        warning('off','stats:jbtest:PTooBig'); % suspend this warning
        warning('off','stats:jbtest:PTooSmall'); % suspend this warning
        [~, p] = jbtest(x);
        warning('on','stats:jbtest:PTooBig'); % resume this warning
        warning('on','stats:jbtest:PTooSmall'); % resume this warning
        
    case 'lbq' % Econometrics Toolbox
        %% Check that an Econometrics Toolbox license is available:
        BF_CheckToolbox('econometrics_toolbox')
        
        % Perform the test
        [~, p] = lbqtest(x);
        
    otherwise
        error('Unknown hypothesis test ''%s''',thetest);
end

end