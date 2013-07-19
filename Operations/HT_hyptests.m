% HT_hyptests
% 
% Outputs the p-value from a statistical hypothesis test applied to the
% time series.
% 
% Tests are implemented as functions in MATLAB's Statistics Toolbox.
% (except Ljung-Box Q-test, which uses the Econometrics Toolbox)
% 
% INPUTS,
% x, the input time series
% thetest, the hypothesis test to perform:
%           (i) sign test ('signtest'),
%           (ii) runs test ('runstest'),
%           (iii) variance test ('vartest'),
%           (iv) Z-test ('ztest'),
%           (v) Wilcoxon signed rank test for a zero median ('signrank'),
%           (vi) Jarque-Bera test of composite normality ('jbtest').
%           (vii) Ljung-Box Q-test for residual autocorrelation ('lbq')

function p = HT_hyptests(x,thetest)
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
        [~, p] = jbtest(x);
        
    case 'lbq' % Econometrics Toolbox
        [~, p] = lbqtest(x);
        
    otherwise
        error('Unknown hypothesis test ''%s''',thetest);
end

end