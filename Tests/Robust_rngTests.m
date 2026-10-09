classdef Robust_rngTests < matlab.unittest.TestCase
    % Tests for the portable random-number generator BF_Random and for the
    % features that use it or were made deterministic (CO_AddNoise, FC_Surprise,
    % BF_TieBreakNoise).

    methods (Test)
        function randomMatchesReferenceGenerator(tc)
            % Values from L'Ecuyer's reference implementation (RngStream.c) for the
            % default state (12345 in all six words).
            u = BF_Random(1000000, 12345 * ones(1, 6));
            tc.verifyEqual(u(1:3), [0.12701112204657714; 0.3185275653967945; 0.30918601558327008], 'AbsTol', 0);
            tc.verifyEqual(u(1000), 0.98607848680213228, 'AbsTol', 0);
            tc.verifyEqual(u(100000), 0.69628910995743587, 'AbsTol', 0);
            tc.verifyEqual(u(1000000), 0.37578835621568801, 'AbsTol', 0);
            % and for the state of all ones
            v = BF_Random(3, ones(1, 6));
            tc.verifyEqual(v, [0.0003395772237870988; 0.55588071598279964; 0.014204660652803588], 'AbsTol', 0);
        end

        function randomIsReproducibleAndLeavesGlobalStreamAlone(tc)
            rng(5); a = rand;
            rng(5); x1 = BF_Random(100, 3, 'normal'); b = rand;
            x2 = BF_Random(100, 3, 'normal');
            tc.verifyEqual(a, b); % global stream not consumed
            tc.verifyEqual(x1, x2);
            tc.verifyNotEqual(x1, BF_Random(100, 4, 'normal'));
        end

        function randomHasRightDistributions(tc)
            u = BF_Random(20000, 1);
            tc.verifyGreaterThan(min(u), 0);
            tc.verifyLessThan(max(u), 1);
            tc.verifyEqual(mean(u), 0.5, 'AbsTol', 0.01);
            z = BF_Random(20001, 2, 'normal');
            tc.verifyEqual(numel(z), 20001);
            tc.verifyEqual(mean(z), 0, 'AbsTol', 0.03);
            tc.verifyEqual(std(z), 1, 'AbsTol', 0.03);
            p = BF_Random(50, 7, 'perm');
            tc.verifyEqual(sort(p), (1:50)');
            tc.verifyNotEqual(p, (1:50)');
            tc.verifyEqual(BF_Random(1, 0, 'perm'), 1);
        end

        function addNoiseIsDeterministicAndStable(tc)
            rng(1); y = zscore(filter(1, [1 -0.8], randn(500, 1)));
            o1 = CO_AddNoise(y, 1, 'even', 10);
            o2 = CO_AddNoise(y, 1, 'even', 10);
            tc.verifyEqual(o1, o2);
            tc.verifyTrue(o1.fitexpa >= 0 && o1.fitexpr2 >= 0 && o1.fitexpr2 <= 1);
            % invariant (to a tiny tolerance) to a 1e-9 perturbation of the series
            rng(2); yp = y + 1e-9 * randn(size(y));
            o3 = CO_AddNoise(yp, 1, 'even', 10);
            tc.verifyEqual(o3.ami_at_10, o1.ami_at_10, 'AbsTol', 0.01);
            tc.verifyEqual(o3.fitexpb, o1.fitexpb, 'AbsTol', 0.05);
        end

        function surpriseIsDeterministicAndBounded(tc)
            rng(3); y = randn(1000, 1);
            o1 = FC_Surprise(y, 'dist', 20, 3, 'quantile', 500, 'default');
            rng(99);
            o2 = FC_Surprise(y, 'dist', 20, 3, 'quantile', 500, 'other');
            tc.verifyEqual(o1, o2);
            tc.verifyGreaterThanOrEqual(o1.effectSize, 0);
            % a constant-surprise setting gives NaN, not a rounding-noise ratio
            ys = repmat([1; 2; 3], 400, 1);
            o = FC_Surprise(ys, 'T1', 20, 3, 'quantile', 500, 'default');
            tc.verifyTrue(isnan(o.effectSize) || o.effectSize < 1e6);
            % long series: subsample is deterministic too
            yl = randn(5000, 1);
            tc.verifyEqual(FC_Surprise(yl, 'dist', 50, 3, 'quantile', 500).mean, ...
                           FC_Surprise(yl, 'dist', 50, 3, 'quantile', 500).mean);
        end

        function tieBreakNoiseIsDeterministic(tc)
            rng(4); yt = round(3 * randn(300, 1));
            a = BF_TieBreakNoise(yt);
            rng(77); b = BF_TieBreakNoise(yt);
            tc.verifyEqual(a, b);
            tc.verifyEqual(numel(unique(a)), numel(a)); % ties are broken
            tc.verifyLessThan(max(abs(a - yt)), 1e-8);
            yu = randn(300, 1); % untied: unchanged
            tc.verifyEqual(BF_TieBreakNoise(yu), yu);
        end
    end
end
