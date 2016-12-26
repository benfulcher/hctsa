function p = HT_DistributionTest(x,theTest,theDistn,numBins)
% HT_DistributionTest   Hypothesis test for distributional fits to a data vector.
%
% Fits a distribution to the data and then performs an appropriate hypothesis
% test to quantify the difference between the two distributions.
%
% We fit Gaussian, Extreme Value, Uniform, Beta, Rayleigh, Exponential, Gamma,
% Log-Normal, and Weibull distributions, using code described for DN_M_kscomp.
%
%---INPUTS:
% x, the input data vector
% theTest, the hypothesis test to perform:
%           (i) 'chi2gof': chi^2 goodness of fit test
%           (ii) 'ks': Kolmogorov-Smirnov test
%           (iii) 'lillie': Lilliefors test
%
% theDistn, the distribution to fit:
%           (i) 'norm' (Normal)
%           (ii) 'ev' (Extreme value)
%           (iii) 'uni' (Uniform)
%           (iv) 'beta' (Beta)
%           (v) 'rayleigh' (Rayleigh)
%           (vi) 'exp' (Exponential)
%           (vii) 'gamma' (Gamma)
%           (viii) 'logn' (Log-normal)
%           (ix) 'wbl' (Weibull)
%
% numBins, the number of bins to use for the chi2 goodness of fit test
%
% All of these functions for hypothesis testing are implemented in Matlab's
% Statistics Toolbox.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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
%% First fit the distribution
% ------------------------------------------------------------------------------
switch theDistn
    case 'norm'
        [a, b] = normfit(x);
    case 'ev'
        a = evfit(x);
    case 'uni'
        [a, b] = unifit(x);
    case 'beta'
        % clumsily scale to the range (0,1)
        x = (x - min(x) + 0.01*std(x)) / (max(x) - min(x) + 0.01*std(x));
        a = betafit(x);        % then fit
    case 'rayleigh'
        if any(x < 0)
            p = NaN; return
        else % valid domain to fit a Rayleigh distribution
            a = raylfit(x);
        end
    case 'exp'
        if any(x < 0)
            p = NaN; return
        else
            a = expfit(x);
        end
    case 'gamma'
        if any(x < 0)
            p = NaN; return
        else
            a = gamfit(x);
        end
    case 'logn'
        if any(x<=0)
            p = NaN; return
        else
            a = lognfit(x);
        end
    case 'wbl'
        if any(x<=0)
            p = NaN; return
        else
            a = wblfit(x);
        end
    otherwise
        error('Unknown distibution ''%s''',theDistn);
end

switch theTest
    case 'chi2gof' % PERFORM A CHI2 GOODNESS OF FIT TEST
        switch theDistn
            case 'norm'
                mycdf = {@normcdf,a,b};
            case 'ev'
                mycdf = {@evcdf,a(1),a(2)};
            case 'uni'
                mycdf = {@unifcdf,a,b};
            case 'beta'
                mycdf = {@betacdf,a(1),a(2)};
            case 'rayleigh'
                mycdf = {@raylcdf,a};
            case 'exp'
                mycdf = {@expcdf,a};
            case 'gamma'
                mycdf = {@gamcdf,a(1),a(2)};
            case 'gp'
%                 mycdf = {@gpcdf,a(1),a(2),theta};
            case 'logn'
                mycdf = {@logncdf,a(1),a(2)};
            case 'wbl'
                mycdf = {@wblcdf,a(1),a(2)};
        end
        warning('off','stats:chi2gof:LowCounts') % temporarily disable this warning
        [~, p] = chi2gof(x,'cdf',mycdf,'nbins',numBins);
        warning('on','stats:chi2gof:LowCounts') % temporarily disable this warning

    case 'ks' % KOLMOGOROV-SMIRNOV TEST
        x1 = sort(unique(round(x*1E6)/1E6)); % unique values to 6 places
        if x1(1) > min(x), x1 = [min(x); x1]; end
        if x1(end) < max(x), x1 = [x1; max(x)]; end
        if size(x1,1) < size(x1,2);
            x1 = x1';
        end
        switch theDistn
            case 'norm'
                mycdf = [x1, normcdf(x1,a,b)];
            case 'ev'
                mycdf = [x1, evcdf(x1,a(1),a(2))];
            case 'uni'
                mycdf = [x1, unifcdf(x1,a,b)];
            case 'beta'
                mycdf = [x1, betacdf(x1,a(1),a(2))];
            case 'rayleigh'
                mycdf = [x1, raylcdf(x1,a)];
            case 'exp'
                mycdf = [x1, expcdf(x1,a)];
            case 'gamma'
                mycdf = [x1, gamcdf(x1,a(1),a(2))];
            case 'logn'
                mycdf = [x1, logncdf(x1,a(1),a(2))];
            case 'wbl'
                mycdf = [x1, wblcdf(x1,a(1),a(2))];
        end

        [~, p] = kstest(x,mycdf);

    case 'lillie' % LILLIEFORS TEST
        % Temporarily suspend low/high tabulated p-value warnings that often occur with this hypothesis test
        warning('off','stats:lillietest:OutOfRangePLow'); warning('off','stats:lillietest:OutOfRangePHigh');
        if any(ismember({'norm','ev'},theDistn))
            [~, p] = lillietest(x,0.05,theDistn);
        elseif strcmp('exp',theDistn)
            if any(x < 0)
                p = NaN; return
            else
                [~, p] = lillietest(x,0.05,theDistn);
            end
        else
           p = NaN;
           fprintf(1,'***RETURNED AN UNEXPECTED NAN FOR LILLIEFORS TEST\n');
        end
        warning('on','stats:lillietest:OutOfRangePLow'); warning('on','stats:lillietest:OutOfRangePHigh');
    otherwise
        error('Unknown test ''%s''',theTest);
end

end
