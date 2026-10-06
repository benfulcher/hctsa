classdef Robust_gpSplitTests < matlab.unittest.TestCase
    % Tests for MF_GP_LocalPrediction's deterministic 'spreadgap' splits and its robust
    % tail summaries over the fits.

    methods (Static)
        function y = series(n, seed)
            % an AR(1)-like series made without MATLAB's random streams
            u = BF_Random(n, seed, 'normal');
            y = zscore(filter(1, [1 -0.7], u));
        end
    end

    methods (Test)
        function spreadgapIgnoresTheSeed(tc)
            y = Robust_gpSplitTests.series(300, 3);
            cf = {'covSum', {'covSEiso', 'covNoise'}};
            a = MF_GP_LocalPrediction(y, cf, 10, 3, 4, 'spreadgap', 0, 2);
            b = MF_GP_LocalPrediction(y, cf, 10, 3, 4, 'spreadgap', 7, 2);
            tc.verifyTrue(isstruct(a));
            tc.verifyEqual(a, b);
        end

        function spreadgapCyclesThroughTheTestSets(tc)
            % 2 windows x 8 splits = 16 fits, more than the nchoosek(5, 2) = 10 possible
            % test sets of a 5-sample window: the list is cycled through
            y = Robust_gpSplitTests.series(200, 4);
            cf = {'covSum', {'covSEiso', 'covNoise'}};
            o = MF_GP_LocalPrediction(y, cf, 3, 2, 2, 'spreadgap', [], 8);
            tc.verifyTrue(isstruct(o));
            tc.verifyTrue(isfinite(o.q90abs_run));
        end

        function tailSummariesAreOrdered(tc)
            y = Robust_gpSplitTests.series(400, 5);
            cf = {'covSum', {'covSEiso', 'covNoise'}};
            o = MF_GP_LocalPrediction(y, cf, 10, 3, 6, 'spreadgap', [], 4);
            tc.verifyLessThanOrEqual(o.meanabs_run, o.q90abs_run + 10 * eps);
            tc.verifyLessThanOrEqual(o.q90abs_run, o.maxabs_run);
            tc.verifyLessThanOrEqual(o.q90abs_std_run, o.maxabs_std_run);
            tc.verifyLessThanOrEqual(o.minabs_std_run, o.low25abs_std_run);
            tc.verifyLessThanOrEqual(o.low25abs_std_run, o.meanabs_std_run);
            tc.verifyLessThanOrEqual(o.meanerrbar, o.high25errbar);
            tc.verifyLessThanOrEqual(o.high25errbar, o.maxerrbar);
            tc.verifyLessThanOrEqual(o.minnlml, o.q90nlml);
            tc.verifyLessThanOrEqual(o.q90nlml, o.maxnlml);
        end

        function tailMeanOfQuarter(tc)
            % with 4 windows and 1 split, the lowest quarter is one window: low25abs_std_run
            % is then the minimum
            y = Robust_gpSplitTests.series(300, 6);
            cf = {'covSum', {'covSEiso', 'covNoise'}};
            o = MF_GP_LocalPrediction(y, cf, 10, 3, 4, 'frombefore');
            tc.verifyEqual(o.low25abs_std_run, o.minabs_std_run);
        end
    end
end
