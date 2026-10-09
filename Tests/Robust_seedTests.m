classdef Robust_seedTests < matlab.unittest.TestCase
    % Tests for the features whose dependence on a random draw was removed (deterministic
    % low-discrepancy subsamples) or reduced (averaging over several draws).

    methods (Static)
        function y = series(n, seed)
            % an AR(1)-like series made without MATLAB's random streams
            u = BF_Random(n, seed, 'normal');
            y = zscore(filter(1, [1 -0.7], u));
        end
    end

    methods (Test)
        function spreadPermIsASpreadPermutation(tc)
            p = BF_SpreadPerm(1000);
            tc.verifyEqual(sort(p), (1:1000)');
            tc.verifyEqual(p, BF_SpreadPerm(1000));
            for k = [10 64 200]
                g = diff([0; sort(p(1:k)); 1000]);
                tc.verifyLessThan(max(g), 3 * 1000 / k); % no big gaps in any prefix
            end
            tc.verifyEqual(BF_SpreadPerm(1), 1);
        end

        function deterministicSubsampleFeaturesIgnoreTheSeed(tc)
            y = Robust_seedTests.series(800, 4);
            a = NL_DelayTime(y, {'ac', 10}, {'ac', 1}, 0);
            b = NL_DelayTime(y, {'ac', 10}, {'ac', 1}, 7);
            tc.verifyEqual(a, b);
            a = SY_SpreadRandomLocal(y, 100, 50, 0);
            b = SY_SpreadRandomLocal(y, 100, 50, 7);
            tc.verifyEqual(a, b);
            a = NL_RQA(y, 1, 3, {'ac', 1}, 0.1, 2, 2, 'full', 0);
            b = NL_RQA(y, 1, 3, {'ac', 1}, 0.1, 2, 2, 'full', 7);
            tc.verifyEqual(a, b);
            a = NL_FractalDimensions(y, 2, 20, 0.2, 1, 5, {'ac', 1}, 8, {1, 3}, 0);
            b = NL_FractalDimensions(y, 2, 20, 0.2, 1, 5, {'ac', 1}, 8, {1, 3}, 7);
            tc.verifyEqual(a, b);
            a = NL_EVTLocalDim(y, 1, 3, 0.95, {'ac', 1}, 50, 5, 'full', 0);
            b = NL_EVTLocalDim(y, 1, 3, 0.95, {'ac', 1}, 50, 5, 'full', 7);
            tc.verifyEqual(a, b);
        end

        function spreadSegmentsAreDeterministic(tc)
            y = Robust_seedTests.series(600, 6);
            a = MF_FitSubsegments(y, 'ar', 2, 'spread', [20, 0.1]);
            b = MF_FitSubsegments(y, 'ar', 2, 'spread', [20, 0.1]);
            tc.verifyEqual(a, b);
            c = MF_CompareTestSets(y, 'ar', 2, 'spread', [20, 0.1], 1, 0);
            d = MF_CompareTestSets(y, 'ar', 2, 'spread', [20, 0.1], 1, 9);
            tc.verifyEqual(c, d);
        end

        function randomizeAveragesRepeats(tc)
            y = Robust_seedTests.series(300, 2);
            % one repeat is the previous single draw; averaging reduces the dependence on the seed
            a1 = EN_Randomize(y, 'statdist', 0, 1); b1 = EN_Randomize(y, 'statdist', 1, 1);
            a20 = EN_Randomize(y, 'statdist', 0, 20); b20 = EN_Randomize(y, 'statdist', 1, 20);
            tc.verifyLessThan(abs(a20.statav5diff - b20.statav5diff), abs(a1.statav5diff - b1.statav5diff));
            tc.verifyEqual(EN_Randomize(y, 'permute', 3), EN_Randomize(y, 'permute', 3)); % reproducible
            % the first checkpoint is the original series in every repeat
            tc.verifyEqual(a20.d1diff >= 0, true);
        end

        function gpHyperparametersAveragesDraws(tc)
            y = Robust_seedTests.series(400, 5);
            cf = {'covSum', {'covSEiso', 'covNoise'}};
            one = MF_GP_Hyperparameters(y, cf, 1, 50, 'random_i', 0, 1);
            two = MF_GP_Hyperparameters(y, cf, 1, 50, 'random_i', 2, 1);
            avg = MF_GP_Hyperparameters(y, cf, 1, 50, 'random_i', 0, 2); % draws with seeds 0 and 2
            tc.verifyEqual(avg.logh1, (one.logh1 + two.logh1) / 2, 'AbsTol', 1e-12);
            tc.verifyEqual(avg.nlml, (one.nlml + two.nlml) / 2, 'AbsTol', 1e-12);
        end

        function gpLocalPredictionSplits(tc)
            y = Robust_seedTests.series(400, 8);
            cf = {'covSum', {'covSEiso', 'covNoise'}};
            a = MF_GP_LocalPrediction(y, cf, 10, 3, 6, 'randomgap', 1, 1);
            b = MF_GP_LocalPrediction(y, cf, 10, 3, 6, 'randomgap', 1, 3);
            tc.verifyTrue(isstruct(a) && isstruct(b));
            tc.verifyEqual(MF_GP_LocalPrediction(y, cf, 10, 3, 6, 'randomgap', 1, 3), b); % reproducible
            % other modes are unaffected by numSplits
            tc.verifyEqual(MF_GP_LocalPrediction(y, cf, 10, 3, 6, 'frombefore', 1, 4), MF_GP_LocalPrediction(y, cf, 10, 3, 6, 'frombefore', 1, 1));
        end

        function c1UsesAtLeast500Centers(tc)
            y = Robust_seedTests.series(1000, 9);
            a = NL_c1(y, 1, [2, 4], 25, 0.1); % 100 requested: raised to 500
            b = NL_c1(y, 1, [2, 4], 25, 500);
            tc.verifyEqual(a, b);
        end
    end
end
