classdef Robust_garchTests < matlab.unittest.TestCase
    % Tests for the GARCH estimator (BF_GARCHFit) and the operations that use it
    % (MF_GARCHfit, MF_GARCHcompare): the likelihood convention, the gradient, convergence,
    % parameter recovery, the test for conditional heteroskedasticity, and the three model types.

    methods(TestClassSetup)
        function runStartup(testCase)
            try
                run("../startup.m")
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'HCTSA failed to startup successfully.')
        end
    end

    methods(Static)
        function y = simulate(n, omega, alpha, beta, seed, gamma, nu)
            % GARCH(1,1) / GJR(1,1) with Gaussian or unit-variance Student-t innovations
            rng(seed);
            if nargin < 6, gamma = 0; end
            if nargin < 7, nu = []; end
            m = n + 500;
            if isempty(nu)
                z = randn(m, 1);
            else
                z = trnd(nu, m, 1) * sqrt((nu - 2) / nu);
            end
            y = zeros(m, 1);
            s2 = omega / (1 - alpha - gamma/2 - beta);
            ePrev = 0;
            for t = 1:m
                s2 = omega + (alpha + gamma*(ePrev < 0)) * ePrev^2 + beta * s2;
                y(t) = sqrt(s2) * z(t);
                ePrev = y(t);
            end
            y = zscore(y(501:end));
        end
    end

    methods(Test)

        function likelihoodConvention(testCase)
            % Gaussian zero-mean GARCH(1,1) log-likelihood with presample variance and
            % squared innovation both equal to the mean of the squared series, as a plain loop
            y = Robust_garchTests.simulate(1500, 0.1, 0.1, 0.8, 3);
            D = BF_GARCHFit('data', y);
            f = BF_GARCHFit('mle', D, 1, 1, false, false);
            w = f.theta(1); b = f.theta(2); a = f.theta(3);
            v0 = mean(y.^2); s2 = v0; e2 = v0; llf = 0;
            for t = 1:length(y)
                s2 = w + b*s2 + a*e2;
                llf = llf - 0.5*(log(2*pi) + log(s2) + y(t)^2/s2);
                e2 = y(t)^2;
            end
            testCase.verifyEqual(f.llf, llf, 'AbsTol', 1e-8);
            testCase.verifyEqual(f.h(1), w + (a + b)*v0, 'RelTol', 1e-12);
        end

        function parameterRecovery(testCase)
            y = Robust_garchTests.simulate(6000, 0.1, 0.1, 0.8, 1);
            D = BF_GARCHFit('data', y);
            f = BF_GARCHFit('mle', D, 1, 1, false, false);
            testCase.verifyEqual(f.theta(3), 0.1, 'AbsTol', 0.04);
            testCase.verifyEqual(f.theta(2), 0.8, 'AbsTol', 0.08);
            testCase.verifyLessThan(f.dec, 1e-8);
            % the fit is a maximum: perturbed parameters do not improve it
            S = f.spec; rng(1);
            for k = 1:20
                th = f.theta .* (1 + 0.02*randn(size(f.theta)));
                if all(S.c - S.A*th > 0)
                    h = BF_GARCHFit('mle', D, 1, 1, false, false, th, false);
                    testCase.verifyEqual(h.llf, f.llf, 'AbsTol', 1e-6);
                end
            end
        end

        function whiteNoiseIsConstantVariance(testCase)
            % Without conditional heteroskedasticity the GARCH coefficient is not identified:
            % the constant-variance model is reported
            for seed = 1:3
                rng(seed); y = randn(1000, 1);
                out = MF_GARCHfit(y, 'none', 1, 1);
                outG = MF_GARCHfit(y, 'none', 1, 1, 'default', 'gjr', 'gaussian');
                outT = MF_GARCHfit(rand(1000, 1), 'none', 1, 1, 'default', 'garch', 't'); % (uniform: tails lighter than Gaussian)
                testCase.verifyEqual(out.hasStructure, 0);
                testCase.verifyEqual([out.ARCH_1, out.GARCH_1, out.persistence, out.stdsigma], [0 0 0 0]);
                testCase.verifyTrue(isnan(out.ARCHerr_1));
                testCase.verifyEqual(outG.leverage, 0);
                testCase.verifyEqual(outT.distDoF, 100);
                testCase.verifyTrue(isfinite(out.uncondVar));
            end
        end

        function structureFoundAndOutputs(testCase)
            y = Robust_garchTests.simulate(1500, 0.1, 0.15, 0.75, 5);
            out = MF_GARCHfit(y, 'none', 1, 1);
            testCase.verifyEqual(out.hasStructure, 1);
            testCase.verifyGreaterThan(out.persistence, 0.4);
            testCase.verifyLessThan(out.persistence, 1);
            testCase.verifyEqual(out.persistence, out.GARCH_1 + out.ARCH_1, 'AbsTol', 1e-12);
            testCase.verifyTrue(isfinite(out.ARCHerr_1));
            outT = MF_GARCHfit(y, 'none', 1, 1, 'default', 'garch', 't');
            outG = MF_GARCHfit(y, 'none', 1, 1, 'default', 'gjr', 'gaussian');
            testCase.verifyGreaterThanOrEqual(outT.LLF, out.LLF - 1e-3);
            testCase.verifyGreaterThanOrEqual(outG.LLF, out.LLF - 1e-6);
            testCase.verifyError(@() MF_GARCHfit(y, 'none', 1, 1, 'default', 'egarch', 'gaussian'), ?MException);
        end

        function leverageAndTailsRecovered(testCase)
            y = Robust_garchTests.simulate(8000, 0.03, 0.03, 0.88, 4, 0.12);
            out = MF_GARCHfit(y, 'none', 1, 1, 'default', 'gjr', 'gaussian');
            testCase.verifyEqual(out.hasStructure, 1);
            testCase.verifyEqual(out.leverage, 0.12, 'AbsTol', 0.06);
            y = Robust_garchTests.simulate(8000, 0.02, 0.08, 0.9, 5, 0, 5);
            out = MF_GARCHfit(y, 'none', 1, 1, 'default', 'garch', 't');
            testCase.verifyEqual(out.distDoF, 5, 'AbsTol', 1.5);
        end

        function reproducibleUnderTinyPerturbation(testCase)
            y = Robust_garchTests.simulate(1000, 0.1, 0.12, 0.8, 7);
            rng(0);
            a = MF_GARCHfit(y, 'none', 1, 1);
            b = MF_GARCHfit(y + 1e-9*randn(size(y)), 'none', 1, 1);
            for f = {'ARCH_1', 'GARCH_1', 'persistence', 'LLF', 'maxsigma'}
                testCase.verifyEqual(b.(f{1}), a.(f{1}), 'AbsTol', 1e-5);
            end
        end

        function compareNestedAndConstantVarianceChoice(testCase)
            y = Robust_garchTests.simulate(1200, 0.1, 0.15, 0.75, 6);
            out = MF_GARCHcompare(y, 'none', 1:2, 1:2);
            testCase.verifyEqual(out.hasStructure, 1);
            testCase.verifyGreaterThanOrEqual(out.maxLLF, out.minLLF); % nested models
            testCase.verifyTrue(ismember(out.bestpBIC, 0:2) && ismember(out.bestqBIC, 0:2));
            rng(3); w = MF_GARCHcompare(randn(1000, 1), 'none', 1:2, 1:2);
            testCase.verifyEqual(w.hasStructure, 0);
            testCase.verifyEqual(w.minK, w.maxK); % constant variance
            testCase.verifyEqual(w.Ks_vary_p, 0);
            testCase.verifyEqual([w.bestpAIC, w.bestqAIC, w.bestpBIC, w.bestqBIC], [0 0 0 0]);
        end

        function pvalueTestsMatchEconometricsToolbox(testCase)
            testCase.assumeTrue(license('test', 'econometrics_toolbox') == 1);
            rng(2); x = randn(300, 1) .* (1 + 0.5*sin((1:300)'/20));
            [~, p1] = archtest(x, 'lags', 1:20);
            testCase.verifyEqual(BF_ARCHTest(x, 20), p1(:)', 'AbsTol', 1e-10);
            [~, q1] = lbqtest(x.^2, 'lags', 1:20);
            testCase.verifyEqual(BF_LjungBoxTest(x.^2, 20), q1(:)', 'AbsTol', 1e-10);
        end
    end
end
