classdef MinLengthTests < matlab.unittest.TestCase
    % Tests of the minimum time-series lengths of the hctsa feature set
    % (FeatureSets/INP_minlength_hctsa.txt, read by TS_AddMinLengths and applied
    % in TS_CalculateFeatureVector):
    %   (1) every bound refers to a master operation or operation of the set;
    %   (2) at its bound, a master operation runs without error on white noise and
    %       an AR(1) process, and one sample below it all its operations are NaN
    %       (without the master operation being evaluated);
    %   (3) an operation with its own bound is NaN one sample below it;
    %   (4) the full feature set gives no errors on short series, with the same
    %       operations as on a long one;
    %   (5) series shorter than 3 samples and constant series give all NaN;
    %   (6) tables without bounds (custom feature sets, older HCTSA files) are
    %       computed at any length, and TS_Init/TS_Compute store and apply bounds.

    properties
        Operations
        MasterOperations
        boundFile = 'INP_minlength_hctsa.txt';
    end

    methods (Static)
        function x = testSeries(whichOne, N)
            % Fixed-seed white noise and AR(1) (phi = 0.6) of length N
            rng(20261103);
            e = randn(N + 100, 1);
            switch whichOne
            case 'noise'
                x = e(1:N);
            case 'ar1'
                x = filter(1, [1, -0.6], e);
                x = x(101:end);
            end
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
            oldState = get(groot, 'defaultFigureVisible');
            set(groot, 'defaultFigureVisible', 'off');
            testCase.addTeardown(@() set(groot, 'defaultFigureVisible', oldState));
        end

        function loadLibrary(testCase)
            ops = TS_ReadInputFile('ops', 'INP_ops_hctsa.txt', false);
            mops = TS_ReadInputFile('mops', 'INP_mops_hctsa.txt', false);
            [ops, mops] = TS_LinkOperationsWithMasters(ops, mops);
            [testCase.Operations, testCase.MasterOperations] = TS_AddMinLengths(ops, mops, 'hctsa');
        end
    end

    methods (Access = private)
        function [fv, cq] = computeMaster(testCase, x, masterInd)
            % Compute the operations of one master operation on x
            mops = testCase.MasterOperations(masterInd, :);
            ops = testCase.Operations(testCase.Operations.MasterID == mops.ID, :);
            [fv, ~, cq] = TS_CalculateFeatureVector(x, false, ops, mops, false, 'fast');
        end

        function x = seriesForBound(~, whichOne, N, isDecimated)
            % A series whose relevant length is N. Master operations on the
            % decimated series are tested on white noise, which decimation leaves
            % unchanged (its autocorrelation falls below 1/e at lag 1).
            if isDecimated
                whichOne = 'noise';
            end
            x = MinLengthTests.testSeries(whichOne, N);
            if isDecimated
                assert(numel(BF_PreProcess(zscore(x), 'decimate_ac1e')) == N, ...
                    'The test series is changed by decimation.');
            end
        end
    end

    methods (Test)
        function test_BoundKeysMatchFeatureSet(testCase)
            fid = fopen(testCase.boundFile);
            dataIn = textscan(fid, '%s%f', 'CommentStyle', '#');
            fclose(fid);
            keys = dataIn{1};
            isOpKey = contains(keys, '.');
            testCase.verifyEmpty(setdiff(keys(~isOpKey), testCase.MasterOperations.Label), ...
                'Master-operation bounds with no matching master operation.');
            testCase.verifyEmpty(setdiff(keys(isOpKey), testCase.Operations.CodeString), ...
                'Operation bounds with no matching operation.');
            testCase.verifyEqual(nnz(testCase.MasterOperations.MinLength) + ...
                nnz(testCase.Operations.MinLength), numel(keys));
            % An operation's own bound exceeds its master operation's:
            withOpBound = find(testCase.Operations.MinLength > 0);
            for i = withOpBound'
                masterBound = testCase.MasterOperations.MinLength(testCase.MasterOperations.ID == testCase.Operations.MasterID(i));
                testCase.verifyGreaterThan(testCase.Operations.MinLength(i), masterBound, ...
                    sprintf('%s: bound not above its master operation''s', testCase.Operations.CodeString{i}));
            end
        end

        function test_MasterOperationsAtTheirBound(testCase)
            mops = testCase.MasterOperations;
            withBound = find(mops.MinLength > 0);
            for k = withBound'
                L = mops.MinLength(k);
                isDecimated = contains(mops.Code{k}, "BF_PreProcess(x_z,'decimate_ac1e')");
                for whichOne = {'noise', 'ar1'}
                    % At the bound: computed without error
                    x = testCase.seriesForBound(whichOne{1}, max(L, 3), isDecimated);
                    [~, cq] = testCase.computeMaster(x, k);
                    testCase.verifyFalse(any(cq == 1), sprintf('%s errors at its bound, N = %u (%s)', ...
                        mops.Code{k}, L, whichOne{1}));
                    % One sample below: all NaN
                    if L > 3
                        x = testCase.seriesForBound(whichOne{1}, L - 1, isDecimated);
                        [fv, cq] = testCase.computeMaster(x, k);
                        testCase.verifyTrue(all(isnan(fv) & cq == 2), sprintf(['%s is not all NaN ' ...
                            'one sample below its bound, N = %u (%s)'], mops.Code{k}, L - 1, whichOne{1}));
                    end
                end
            end
        end

        function test_OperationsBelowTheirBound(testCase)
            ops = testCase.Operations;
            withBound = find(ops.MinLength > 0);
            for i = withBound'
                L = ops.MinLength(i);
                mi = find(testCase.MasterOperations.ID == ops.MasterID(i));
                isDecimated = contains(testCase.MasterOperations.Code{mi}, "BF_PreProcess(x_z,'decimate_ac1e')");
                x = testCase.seriesForBound('ar1', L - 1, isDecimated);
                [fv, ~, cq] = TS_CalculateFeatureVector(x, false, ops(i, :), testCase.MasterOperations(mi, :), false, 'fast');
                testCase.verifyTrue(isnan(fv) && cq == 2, sprintf('%s is not NaN at N = %u', ops.CodeString{i}, L - 1));
                x = testCase.seriesForBound('ar1', L, isDecimated);
                [~, ~, cq] = TS_CalculateFeatureVector(x, false, ops(i, :), testCase.MasterOperations(mi, :), false, 'fast');
                testCase.verifyNotEqual(cq, 1, sprintf('%s errors at its bound, N = %u', ops.CodeString{i}, L));
            end
        end

        function test_FullFeatureSetOnShortSeries(testCase)
            % No errors on short series, and every operation is returned
            for N = [3, 8, 15]
                for whichOne = {'noise', 'ar1'}
                    x = MinLengthTests.testSeries(whichOne{1}, N);
                    [fv, ~, cq, masterErrors] = TS_CalculateFeatureVector(x, false, ...
                        testCase.Operations, testCase.MasterOperations, false, 'fast');
                    testCase.verifySize(fv, [height(testCase.Operations), 1]);
                    testCase.verifyEmpty(masterErrors, sprintf('Errors at N = %u (%s):\n%s', N, whichOne{1}, ...
                        strjoin(masterErrors(:, 1), newline)));
                    testCase.verifyFalse(any(cq == 1));
                end
            end
        end

        function test_GlobalFloor(testCase)
            ops = testCase.Operations;
            mops = testCase.MasterOperations;
            for x = {[0.3; -1.2], 2.5 * ones(50, 1)}
                [fv, ~, cq] = TS_CalculateFeatureVector(x{1}, false, ops, mops, false, 'fast');
                testCase.verifyTrue(all(isnan(fv)) && all(cq == 2));
            end
            % The floor applies without bounds too:
            ops(:, 'MinLength') = [];
            mops(:, 'MinLength') = [];
            [fv, ~, cq] = TS_CalculateFeatureVector([1; 2], false, ops(1:20, :), mops, false, 'fast');
            testCase.verifyTrue(all(isnan(fv)) && all(cq == 2));
        end

        function test_NoBoundsComputesAtAnyLength(testCase)
            % Without MinLength columns a master operation is evaluated below its bound
            mops = testCase.MasterOperations;
            k = find(mops.MinLength > 4 & ~contains(mops.Code, 'decimate_ac1e'), 1);
            x = MinLengthTests.testSeries('noise', mops.MinLength(k) - 1);
            ops = testCase.Operations(testCase.Operations.MasterID == mops.ID(k), :);
            [~, ~, cqBound] = TS_CalculateFeatureVector(x, false, ops, mops(k, :), false, 'fast');
            ops(:, 'MinLength') = [];
            mopsNoBound = mops(k, :);
            mopsNoBound(:, 'MinLength') = [];
            [fvFree, ~, cqFree] = TS_CalculateFeatureVector(x, false, ops, mopsNoBound, false, 'fast');
            testCase.verifyTrue(all(cqBound == 2));
            testCase.verifyFalse(all(isnan(fvFree) & cqFree == 2), ...
                sprintf('%s not evaluated without a bound', mops.Code{k}));
        end

        function test_InitAndComputeApplyBounds(testCase)
            % TS_Init stores the bounds; TS_Compute applies them in each mode
            tempDir = tempname;
            mkdir(tempDir);
            testCase.addTeardown(@() rmdir(tempDir, 's'));
            mopsFile = fullfile(tempDir, 'INP_mops_test.txt');
            opsFile = fullfile(tempDir, 'INP_ops_test.txt');
            boundsFile = fullfile(tempDir, 'INP_minlength_test.txt');
            writelines(["CO_trev(x_z,10)	CO_trev_10", "DN_Moments(x_z,3)	DN_Moments_3"], mopsFile);
            writelines(["CO_trev_10.raw	trev_10_raw	test", "CO_trev_10.num	trev_10_num	test", ...
                        "DN_Moments_3	DN_Moments_3	test"], opsFile);
            writelines(["# test bounds", "CO_trev_10	20", "DN_Moments_3	5"], boundsFile);
            timeSeriesData = {MinLengthTests.testSeries('noise', 19), MinLengthTests.testSeries('noise', 30), ...
                              MinLengthTests.testSeries('noise', 4)};
            labels = {'N19', 'N30', 'N4'};
            keywords = {'test', 'test', 'test'};
            inpFile = fullfile(tempDir, 'INP_test.mat');
            save(inpFile, 'timeSeriesData', 'labels', 'keywords');
            for parallelHow = {'none', 'operations'}
                hctsaFile = fullfile(tempDir, sprintf('HCTSA_%s.mat', parallelHow{1}));
                TS_Init(inpFile, {mopsFile, opsFile, boundsFile}, false, hctsaFile);
                loaded = load(hctsaFile, 'Operations', 'MasterOperations');
                testCase.verifyEqual(sort(loaded.MasterOperations.MinLength)', [5, 20]);
                TS_Compute(parallelHow{1}, [], [], 'missing', hctsaFile, 'fast');
                loaded = load(hctsaFile, 'TS_Quality', 'Operations');
                isTrev = startsWith(loaded.Operations.CodeString, 'CO_trev');
                testCase.verifyTrue(all(loaded.TS_Quality(1, isTrev) == 2)); % N = 19 < 20
                testCase.verifyTrue(all(loaded.TS_Quality(2, :) == 0)); % N = 30
                testCase.verifyTrue(all(loaded.TS_Quality(3, :) == 2)); % N = 4 < 5
            end
        end
    end
end
