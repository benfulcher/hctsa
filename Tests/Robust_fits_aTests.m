classdef Robust_fits_aTests < matlab.unittest.TestCase
    % Unit tests for features whose curve fits were made robust: global,
    % deterministic fits (no starting points or optimizer tolerances), bounded
    % goodness-of-fit statistics, and insensitivity to tiny perturbations.

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

    methods(Test)

        %-------------------------------------------------------------
        % BF_TheilSen
        %-------------------------------------------------------------
        function test_BF_TheilSen_ExactLineAndOutliers(testCase)
            x = (1:30)';
            y = 2*x + 3;
            testCase.verifyEqual(BF_TheilSen(x, y), [2, 3], 'AbsTol', 1e-12);
            % a few gross outliers barely move the fit (least squares would):
            y([4 17 25]) = [2000 -1500 3000];
            p = BF_TheilSen(x, y);
            testCase.verifyEqual(p, [2, 3], 'AbsTol', 1e-9);
            pols = polyfit(x, y, 1);
            testCase.verifyGreaterThan(abs(pols(1) - 2), 0.5);
        end

        %-------------------------------------------------------------
        % BF_ExpFit
        %-------------------------------------------------------------
        function test_BF_ExpFit_RecoversParameters(testCase)
            x = linspace(0, 4, 60)';
            f = BF_ExpFit(x, 1.5*exp(-0.9*x) + 0.4, true);
            testCase.verifyEqual([f.a f.b f.c], [1.5 -0.9 0.4], 'AbsTol', 1e-6);
            testCase.verifyEqual(f.r2, 1, 'AbsTol', 1e-9);
            f = BF_ExpFit(x, 2*exp(0.5*x), false);
            testCase.verifyEqual([f.a f.b f.c], [2 0.5 0], 'AbsTol', 1e-6);
        end

        function test_BF_ExpFit_BoundedAndDeterministic(testCase)
            rng(3);
            x = (1:10)';
            for rep = 1:20
                y = randn(10, 1); % pure noise: no exponential structure
                f1 = BF_ExpFit(x, y, true);
                f2 = BF_ExpFit(x, y, true);
                testCase.verifyEqual(f1, f2);
                testCase.verifyGreaterThanOrEqual(f1.r2, 0);
                testCase.verifyLessThanOrEqual(f1.r2, 1);
                testCase.verifyLessThanOrEqual(abs(f1.b), 20/range(x) + 1e-12);
                g = BF_ExpFit(x, y, false);
                testCase.verifyGreaterThanOrEqual(g.r2, 0);
            end
            % a nearly-straight line: r2 stays high even though a and c are not determined
            f = BF_ExpFit(x, 3 + 0.2*x, true);
            testCase.verifyGreaterThan(f.r2, 1 - 1e-9);
            % constant data: NaN
            testCase.verifyTrue(isnan(BF_ExpFit(x, ones(10, 1), true).r2));
        end

        %-------------------------------------------------------------
        % SC_FluctAnal
        %-------------------------------------------------------------
        function test_SC_FluctAnal_Scaling(testCase)
            rng(5);
            wn = randn(2000, 1);
            out = SC_FluctAnal(zscore(wn), 2, 'dfa', 50, 1, [], true);
            testCase.verifyEqual(out.alpha, 0.5, 'AbsTol', 0.1); % white noise
            out = SC_FluctAnal(zscore(cumsum(wn)), 2, 'dfa', 50, 1, [], true);
            testCase.verifyEqual(out.alpha, 1.5, 'AbsTol', 0.15); % random walk
        end

        function test_SC_FluctAnal_BoundedDeterministicStable(testCase)
            rng(6);
            x = zscore(filter(1, [1 -0.8], randn(1000, 1)));
            xp = zscore(x + 1e-9*randn(1000, 1));
            wtfs = {'dfa', 'rsrange', 'std', 'endptdiff'};
            for k = 1:numel(wtfs)
                a = SC_FluctAnal(x, 2, wtfs{k}, 50, 1, [], true);
                b = SC_FluctAnal(x, 2, wtfs{k}, 50, 1, [], true);
                c = SC_FluctAnal(xp, 2, wtfs{k}, 50, 1, [], true);
                testCase.verifyEqual(a, b);
                testCase.verifyGreaterThanOrEqual(a.splitgain, 0);
                testCase.verifyLessThanOrEqual(a.splitgain, 1);
                testCase.verifyEqual(a.alphadiff, a.r1_alpha - a.r2_alpha);
                testCase.verifyFalse(isfield(a, 'alpharat') || isfield(a, 'ratsplitminerr'));
                testCase.verifyEqual(c.alpha, a.alpha, 'AbsTol', 1e-6);
                testCase.verifyEqual(c.prop_r1, a.prop_r1, 'AbsTol', 0.1);
            end
        end

        %-------------------------------------------------------------
        % DN_OutlierInclude
        %-------------------------------------------------------------
        function test_DN_OutlierInclude_FitsBoundedAndStable(testCase)
            rng(7);
            x = zscore(randn(1500, 1));
            xp = zscore(x + 1e-9*randn(1500, 1));
            modes = {'abs', 'pos', 'neg'};
            for k = 1:3
                a = DN_OutlierInclude(x, modes{k}, 0.01);
                c = DN_OutlierInclude(xp, modes{k}, 0.01);
                for f = {'mfexpr2', 'nfexpr2', 'stdrfexpr2', 'nflr2', 'stdrflr2'}
                    testCase.verifyGreaterThanOrEqual(a.(f{1}), 0);
                    testCase.verifyLessThanOrEqual(a.(f{1}), 1);
                end
                testCase.verifyEqual(c.mfexpb, a.mfexpb, 'RelTol', 1e-3, 'AbsTol', 1e-3);
                testCase.verifyEqual(c.nfexpr2, a.nfexpr2, 'AbsTol', 1e-4);
                testCase.verifyEqual(c.nfexprmse, a.nfexprmse, 'AbsTol', 1e-3);
            end
        end

        %-------------------------------------------------------------
        % FC_LoopLocalSimple
        %-------------------------------------------------------------
        function test_FC_LoopLocalSimple_ExpFitBoundedAndStable(testCase)
            rng(8);
            x = zscore(filter(1, [1 -0.7], randn(800, 1)));
            xp = zscore(x + 1e-9*randn(800, 1));
            a = FC_LoopLocalSimple(x, 'mean');
            c = FC_LoopLocalSimple(xp, 'mean');
            testCase.verifyGreaterThanOrEqual(a.sws_fexp_r2, 0);
            testCase.verifyLessThanOrEqual(a.sws_fexp_r2, 1);
            testCase.verifyEqual(c.sws_fexp_r2, a.sws_fexp_r2, 'AbsTol', 1e-5);
            testCase.verifyEqual(c.sws_fexp_rmse, a.sws_fexp_rmse, 'AbsTol', 1e-6);
            testCase.verifyFalse(isfield(a, 'sws_fexp_a') || isfield(a, 'sws_fexp_c'));
        end

        %-------------------------------------------------------------
        % SB_TransitionPAlphabet
        %-------------------------------------------------------------
        function test_SB_TransitionPAlphabet_FitsBoundedAndStable(testCase)
            rng(9);
            x = zscore(filter(1, [1 -0.5], randn(2000, 1)));
            xp = zscore(x + 1e-9*randn(2000, 1));
            a = SB_TransitionPAlphabet(x, 2:20, 1);
            c = SB_TransitionPAlphabet(xp, 2:20, 1);
            for f = {'meandiag', 'maxdiag', 'tr', 'trcov', 'stdeig'}
                r2 = a.([f{1} 'fexp_r2']);
                testCase.verifyGreaterThanOrEqual(r2, 0);
                testCase.verifyLessThanOrEqual(r2, 1);
                testCase.verifyEqual(c.([f{1} 'fexp_b']), a.([f{1} 'fexp_b']), 'AbsTol', 1e-6);
                testCase.verifyEqual(c.([f{1} 'fexp_r2']), r2, 'AbsTol', 1e-6);
            end
        end

    end
end
