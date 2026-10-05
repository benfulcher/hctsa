classdef Robust_fits_bTests < matlab.unittest.TestCase
    % Tests of the deterministic fits and runs-test statistic: BF_RunsZ, BF_FitSinusoids,
    % BF_FitDensityCurve, BF_GaussMix2, and the functions that use them (DN_SimpleFit,
    % SP_SinusoidFit, NW_VisibilityGraph, PP_Compare, PH_Walker, HT_IndependenceTests).

    methods (Test)

        function test_BF_RunsZ_FormulaAndDirection(testCase)
            % 4 values above the median and 4 at or below: alternating gives the maximum number
            % of runs (8), sorted gives the minimum (2)
            alt = [1 -1 1 -1 1 -1 1 -1]';
            n1 = 4; n2 = 4; n = 8;
            mu = 1 + 2*n1*n2/n;
            sigma = sqrt(2*n1*n2*(2*n1*n2 - n)/(n^2*(n - 1)));
            testCase.verifyEqual(BF_RunsZ(alt), (8 - mu)/sigma, 'AbsTol', 1e-12);
            testCase.verifyEqual(BF_RunsZ(sort(alt)), (2 - mu)/sigma, 'AbsTol', 1e-12);
            testCase.verifyGreaterThan(BF_RunsZ(alt), 0);
            testCase.verifyLessThan(BF_RunsZ(sort(alt)), 0);
            testCase.verifyTrue(isnan(BF_RunsZ(ones(10, 1)))); % constant: no test possible
            % strongly structured series: finite, large in magnitude (a p-value would underflow)
            z = BF_RunsZ(cumsum(randn(5000, 1)));
            testCase.verifyTrue(isfinite(z) && z < -20);
        end

        function test_BF_RunsZ_ApproximatelyStandardNormal(testCase)
            rng(1);
            z = zeros(500, 1);
            for i = 1:500, z(i) = BF_RunsZ(randn(200, 1)); end
            testCase.verifyEqual(mean(z), 0, 'AbsTol', 0.15);
            testCase.verifyEqual(std(z), 1, 'AbsTol', 0.15);
        end

        function test_BF_FitSinusoids_RecoversFrequenciesAndIsStable(testCase)
            rng(2);
            N = 600; t = (1:N)';
            y = 2*sin(2*pi*0.031*t + 1) + 1.2*sin(2*pi*0.12*t + 0.3) + 0.8*sin(2*pi*0.33*t) + 0.1*randn(N, 1);
            [~, f] = BF_FitSinusoids(y, 3);
            testCase.verifyEqual(f, [0.031; 0.12; 0.33], 'AbsTol', 5e-4);
            % deterministic, and insensitive to a tiny perturbation
            [yfit1, f1] = BF_FitSinusoids(y, 3);
            [yfit2, f2] = BF_FitSinusoids(y + 1e-9*std(y)*randn(N, 1), 3);
            testCase.verifyEqual(f1, f2, 'AbsTol', 1e-6);
            testCase.verifyEqual(yfit1, yfit2, 'AbsTol', 1e-6);
            % frequencies stay in the allowed range even for a trend
            [~, ftr] = BF_FitSinusoids((1:N)' + randn(N, 1), 1);
            testCase.verifyGreaterThanOrEqual(ftr, 1/(2*N) - 1e-12);
        end

        function test_SP_SinusoidFit_OutputsAndNaN(testCase)
            rng(3);
            y = zscore(sin(2*pi*0.05*(1:300)') + 0.5*randn(300, 1));
            out = SP_SinusoidFit(y, 'sin1');
            testCase.verifyGreaterThan(out.r2, 0.6);
            testCase.verifyTrue(all(isfield(out, {'r2', 'adjr2', 'rmse', 'resAC1', 'resAC2', 'resrunsz'})));
            % a perfect fit leaves only rounding noise in the residuals: not meaningful
            clean = SP_SinusoidFit(sin(2*pi*0.05*(1:300)'), 'sin1');
            testCase.verifyTrue(isnan(clean.resrunsz) && isnan(clean.resAC1));
            % too short for the number of parameters
            testCase.verifyTrue(isnan(SP_SinusoidFit(randn(8, 1), 'sin3')));
        end

        function test_BF_FitDensityCurve_RecoversNoiseFreeCurves(testCase)
            x = linspace(-3, 4, 40)';
            g = @(a, m, s) a*exp(-(x - m).^2/(2*s^2));
            p = g(0.4, 0.5, 0.8);
            testCase.verifyEqual(BF_FitDensityCurve(x, p, 'gauss'), p, 'AbsTol', 1e-6);
            p2 = g(0.3, -1, 0.5) + g(0.25, 2, 0.7);
            testCase.verifyEqual(BF_FitDensityCurve(x, p2, 'gauss2'), p2, 'AbsTol', 1e-5);
            p3 = 0.7*exp(-0.8*x);
            testCase.verifyEqual(BF_FitDensityCurve(x, p3, 'exp'), p3, 'AbsTol', 1e-6);
            xp = linspace(1, 20, 30)';
            p4 = 2*xp.^-1.5;
            testCase.verifyEqual(BF_FitDensityCurve(xp, p4, 'power'), p4, 'AbsTol', 1e-6);
        end

        function test_DN_SimpleFit_DeterministicAndPerturbationStable(testCase)
            rng(4);
            x = zscore([randn(600, 1); 3 + 0.5*randn(400, 1)]);
            xp = x + 1e-9*randn(size(x));
            for model = {'gauss1', 'gauss2', 'exp1'}
                a = DN_SimpleFit(x, model{1}, 'sqrt');
                b = DN_SimpleFit(x, model{1}, 'sqrt');
                c = DN_SimpleFit(xp, model{1}, 'sqrt');
                testCase.verifyEqual(a, b);
                testCase.verifyEqual(a.r2, c.r2, 'AbsTol', 1e-4);
                testCase.verifyEqual(a.rmse, c.rmse, 'AbsTol', 1e-4);
            end
            % two clear peaks are described well by two Gaussians and poorly by one
            testCase.verifyGreaterThan(DN_SimpleFit(x, 'gauss2', 'sqrt').r2, 0.9);
            testCase.verifyLessThan(DN_SimpleFit(x, 'gauss1', 'sqrt').r2, DN_SimpleFit(x, 'gauss2', 'sqrt').r2);
            % power law: positive data only
            testCase.verifyTrue(isnan(DN_SimpleFit(randn(500, 1), 'power1', 'sqrt')));
            pw = DN_SimpleFit(exp(randn(1000, 1)), 'power1', 'sqrt');
            testCase.verifyTrue(isfield(pw, 'resrunsz') && isfinite(pw.r2));
        end

        function test_NW_VisibilityGraph_CollinearPointsAndFewBins(testCase)
            % On a straight line only neighboring points are linked (every point in between
            % blocks the view); rounding in the slopes must not create extra links
            y = 0.1*(1:300)' + 5;
            out = NW_VisibilityGraph(y, 'norm');
            testCase.verifyEqual(out.meank, 2*299/300, 'AbsTol', 1e-12);
            % only two distinct degrees (1 and 2): no distribution fit is meaningful
            testCase.verifyTrue(isnan(out.dgaussk_r2) && isnan(out.dexpk_r2) && isnan(out.dpowerk_resrunsz));
            % a series with a spread of degrees gives finite fits, unchanged by a tiny perturbation
            rng(5);
            y = randn(800, 1);
            a = NW_VisibilityGraph(y, 'norm');
            b = NW_VisibilityGraph(y + 1e-9*randn(800, 1), 'norm');
            for f = {'dgaussk_r2', 'dexpk_r2', 'dpowerk_r2', 'dexpk_resAC1'}
                testCase.verifyTrue(isfinite(a.(f{1})));
                testCase.verifyEqual(a.(f{1}), b.(f{1}), 'AbsTol', 1e-6);
            end
        end

        function test_RunsStatisticInUsers(testCase)
            rng(6);
            x = cumsum(randn(500, 1));
            testCase.verifyEqual(HT_IndependenceTests(x, 'runsz'), BF_RunsZ(x));
            testCase.verifyLessThan(HT_IndependenceTests(x, 'runsz'), -10);
            w = PH_Walker(zscore(x), 'prop', 0.5);
            testCase.verifyTrue(isfield(w, 'res_runsz') && isfinite(w.res_runsz));
            out = PP_Compare(x, 'sin1');
            testCase.verifyTrue(isstruct(out) && isfinite(out.gauss1_kd_resrunsz));
        end

    end
end
