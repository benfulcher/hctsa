classdef RobustnessTests < matlab.unittest.TestCase
    % Whole-library robustness tests: run every registered hctsa feature on a
    % small set of adversarial synthetic time series and check that
    %   (1) no operation throws an error (data-dependent failures must return
    %       NaN, per the NaN-vs-error convention), and
    %   (2) recomputing a series gives bit-identical results (every source of
    %       randomness in an operation must be seeded).
    %
    % These complement OperationsUnitTests.m (known answers for individual
    % functions) and BasicPipelineTests.m (the TS_* pipeline end-to-end): an
    % audit in 2026-09 found a pipeline-stalling loop, an irreproducible
    % estimator and several error()s on legitimate data that 68 passing tests
    % had not caught, but that this style of check caught immediately.
    %
    % The series are chosen to hit the failure modes that ordinary test data
    % misses: heavily repeated values (tie handling), near-constant positive
    % data (scale-dependent loops), a random walk (no ACF zero-crossing), a
    % sparse/spiky series (degenerate distributions) and a short series.
    %
    % Each full-library evaluation takes ~30-60 s, so the whole class runs in
    % a few minutes. See also hctsa_smoke_report.m for a broader, non-test
    % diagnostic (constant features, slowest operations) over the same series.

    properties (Constant)
        seriesNames = {'quantized','nearConstant','randomWalk','spikes','short'};
    end

    properties
        Operations
        MasterOperations
    end

    methods (Static)
        function [ts, name] = adversarialSeries(whichOne)
            % Fixed-seed synthetic series designed to expose fragile operations.
            rng(20260920);
            switch whichOne
            case 'quantized' % many exactly-repeated values
                ts = round(2 * randn(1000, 1));
            case 'nearConstant' % positive, raw (un-z-scored) with a tiny coefficient of variation
                ts = 1 + 1e-8 * randn(1000, 1);
            case 'randomWalk' % nonstationary; no autocorrelation zero-crossing
                ts = cumsum(randn(1000, 1));
            case 'spikes' % mostly near-zero with rare large events
                ts = 0.01 * randn(1000, 1);
                ts(randperm(1000, 20)) = 5 + randn(20, 1);
            case 'short'
                ts = randn(100, 1);
            case 'ar1' % a well-behaved continuous series, for the determinism check
                ts = zeros(1000, 1); e = randn(1000, 1);
                for t = 2:1000, ts(t) = 0.8 * ts(t - 1) + e(t); end
            otherwise
                error('Unknown series ''%s''', whichOne);
            end
            name = whichOne;
        end
    end

    methods (TestClassSetup)
        function runStartup(testCase)
            try
                run("../startup.m")
                pass = true;
            catch
                pass = false;
            end
            testCase.fatalAssertTrue(pass, 'HCTSA failed to startup successfully.')
        end

        function loadLibrary(testCase)
            ops = TS_ReadInputFile('ops', 'INP_ops_hctsa.txt', false);
            mops = TS_ReadInputFile('mops', 'INP_mops_hctsa.txt', false);
            [testCase.Operations, testCase.MasterOperations] = TS_LinkOperationsWithMasters(ops, mops);
        end
    end

    methods (Access = private)
        function [fv, cq] = computeAll(testCase, ts)
            [fv, ~, cq] = TS_CalculateFeatureVector(ts, false, testCase.Operations, ...
                            testCase.MasterOperations, true, 'fast');
        end

        function describeErrors(testCase, cq, seriesName)
            erroredMasters = unique(testCase.Operations.MasterID(cq == 1));
            codes = cell(numel(erroredMasters), 1);
            for i = 1:numel(erroredMasters)
                codes{i} = testCase.MasterOperations.Code{testCase.MasterOperations.ID == erroredMasters(i)};
            end
            testCase.verifyEmpty(erroredMasters, sprintf(['%u operation(s) threw an error on the ''%s'' series ' ...
                '(data-dependent failures should return NaN instead):\n  %s'], ...
                numel(erroredMasters), seriesName, strjoin(codes, newline + "  ")));
        end
    end

    methods (Test)
        function test_NoErrorsOnAdversarialSeries(testCase)
            for i = 1:numel(testCase.seriesNames)
                [ts, name] = RobustnessTests.adversarialSeries(testCase.seriesNames{i});
                [~, cq] = testCase.computeAll(ts);
                testCase.describeErrors(cq, name);
            end
        end

        function test_DeterministicOnQuantizedSeries(testCase)
            % Repeated values are where tie-breaking noise (and any other
            % unseeded randomness) shows up as run-to-run differences.
            testCase.checkDeterminism('quantized');
        end

        function test_DeterministicOnContinuousSeries(testCase)
            testCase.checkDeterminism('ar1');
        end

        function test_EachOperationIndependentOfGlobalRngState(testCase)
            % The two tests above recompute the whole library in file order,
            % which cannot catch an operation that draws from the global random
            % stream without seeding it when an EARLIER operation happens to
            % reset that stream (NL_EmbedCluster passed them for exactly that
            % reason). Here every master operation is evaluated twice from two
            % different global RNG states, in isolation, and must agree.
            [ts, ~] = RobustnessTests.adversarialSeries('ar1');
            x = ts; x_z = zscore(x);
            mops = testCase.MasterOperations;
            bad = {};
            for i = 1:height(mops)
                fn = str2func(['@(x,x_z) ', mops.Code{i}]);
                try
                    rng(1); o1 = fn(x, x_z);
                    rng(2); o2 = fn(x, x_z);
                catch
                    continue % errors are the no-error test's business
                end
                if ~isequaln(o1, o2)
                    bad{end+1} = mops.Code{i}; %#ok<AGROW>
                end
            end
            testCase.verifyEmpty(bad, sprintf(['%u master operation(s) gave different outputs from ' ...
                'different global RNG states (unseeded randomness):\n  %s'], numel(bad), strjoin(bad, newline + "  ")));
        end
    end

    methods (Access = private)
        function checkDeterminism(testCase, whichSeries)
            [ts, name] = RobustnessTests.adversarialSeries(whichSeries);
            [fv1, cq1] = testCase.computeAll(ts);
            [fv2, cq2] = testCase.computeAll(ts);
            % Quality codes must agree, and values must agree wherever the
            % result was a real number (special values are coded in cq):
            differs = (cq1 ~= cq2) | (cq1 == 0 & cq2 == 0 & fv1 ~= fv2);
            badMasters = unique(testCase.Operations.MasterID(differs));
            codes = cell(numel(badMasters), 1);
            for i = 1:numel(badMasters)
                codes{i} = testCase.MasterOperations.Code{testCase.MasterOperations.ID == badMasters(i)};
            end
            testCase.verifyEmpty(badMasters, sprintf(['%u feature(s) across %u master operation(s) changed ' ...
                'between two identical computations of the ''%s'' series (unseeded randomness?):\n  %s'], ...
                sum(differs), numel(badMasters), name, strjoin(codes, newline + "  ")));
        end
    end
end
