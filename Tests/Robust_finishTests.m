classdef Robust_finishTests < matlab.unittest.TestCase
    % Tests for the features that draw their random numbers from the portable
    % generator BF_Random (so they do not depend on MATLAB's random stream), and for
    % the explicit histogram/kernel-density call sites.

    methods (Static)
        function y = series(n, seed)
            % an AR(1)-like series made without MATLAB's random streams
            u = BF_Random(n, seed, 'normal');
            y = zscore(filter(1, [1 -0.7], u));
        end
    end

    methods (Test)
        function fastGeneratorMatchesScalarLoop(tc)
            % long streams (jump-ahead blocks) equal the one-draw-at-a-time recurrence
            n = 7000;
            s = 12345 * ones(1, 6); m1 = 4294967087; m2 = 4294944443;
            u = zeros(n, 1);
            for i = 1:n
                p1 = mod(1403580 * s(2) - 810728 * s(1), m1); s(1:3) = [s(2), s(3), p1];
                p2 = mod(527612 * s(6) - 1370589 * s(4), m2); s(4:6) = [s(5), s(6), p2];
                u(i) = (p1 - p2 + m1 * (p1 <= p2)) * 2.328306549295727688e-10;
            end
            tc.verifyEqual(BF_Random(n, 12345 * ones(1, 6)), u, 'AbsTol', 0);
            % permutation = ranks of the uniforms
            [~, ix] = sort(BF_Random(4000, 5));
            tc.verifyEqual(BF_Random(4000, 5, 'perm'), ix);
        end

        function seedValues(tc)
            tc.verifyEqual(BF_RandomSeed([]), 0);
            tc.verifyEqual(BF_RandomSeed('default'), 0);
            tc.verifyEqual(BF_RandomSeed(7), 7);
            rng(1); a = BF_RandomSeed('none'); rng(1); b = BF_RandomSeed('none');
            tc.verifyEqual(a, b); % drawn from the global stream
        end

        function surrogatesAreReproducibleAndPortable(tc)
            y = Robust_finishTests.series(600, 3);
            for how = {'RP', 'AAFT', 'TFT', 'RandPerm'}
                rng(1); a = SD_MakeSurrogates(y, how{1}, 4, [], 'default');
                rng(2); b = SD_MakeSurrogates(y, how{1}, 4, [], 'default');
                tc.verifyEqual(a, b, how{1}); % independent of MATLAB's stream
                c = SD_MakeSurrogates(y, how{1}, 4, [], 5);
                tc.verifyNotEqual(a, c, how{1});
                tc.verifyNotEqual(a(:, 1), a(:, 2), how{1});
            end
            % properties of the surrogates
            rp = SD_MakeSurrogates(y, 'RP', 3, [], 0);
            tc.verifyEqual(abs(fft(rp(:, 2))), abs(fft(y)), 'AbsTol', 1e-8); % spectrum kept
            aa = SD_MakeSurrogates(y, 'AAFT', 3, [], 0);
            tc.verifyEqual(sort(aa(:, 1)), sort(y), 'AbsTol', 1e-12); % amplitude distribution kept
            rpm = SD_MakeSurrogates(y, 'RandPerm', 3, [], 0);
            tc.verifyEqual(sort(rpm(:, 3)), sort(y));
            % odd length
            ro = SD_MakeSurrogates(y(1:599), 'RP', 2, [], 0);
            tc.verifyEqual(abs(fft(ro(:, 1))), abs(fft(y(1:599))), 'AbsTol', 1e-8);
        end

        function removePointsRandomIsReproducible(tc)
            y = Robust_finishTests.series(500, 4);
            rng(1); a = BF_RemovePoints(y, 'random', 0.2, 'remove', 'default');
            rng(2); b = BF_RemovePoints(y, 'random', 0.2, 'remove', 'default');
            tc.verifyEqual(a, b);
            tc.verifyEqual(numel(a), round(500 * 0.8));
            tc.verifyNotEqual(a, BF_RemovePoints(y, 'random', 0.2, 'remove', 3));
        end

        function randomSubsamplingFeaturesIgnoreMatlabStream(tc)
            y = Robust_finishTests.series(1500, 6);
            calls = {@() SY_SpreadRandomLocal(y, 100, 20, 'default'), ...
                     @() SY_LocalGlobal(y, 'randcg', 50), ...
                     @() MF_FitSubsegments(y, 'ar', 2, 'rand', [10, 0.1], 'default'), ...
                     @() EN_Randomize(y(1:300), 'dyndist', 'default'), ...
                     @() EN_Randomize(y(1:300), 'permute', 'default'), ...
                     @() NL_RQA(y, 2, 3, 10, 0.1, 2, 2, 'full', 'default'), ...
                     @() NL_RecurrenceTimes(y, 2, 3, 10, 0.1, 3, 'full', 'default'), ...
                     @() NL_EVTLocalDim(y, 2, 3, 0.95, 10, 100, 3, 'full', 'default')};
            for i = 1:numel(calls)
                rng(1); a = calls{i}();
                rng(2); b = calls{i}();
                tc.verifyEqual(a, b, sprintf('call %d', i));
            end
            % the seed does matter
            a = SY_SpreadRandomLocal(y, 100, 20, 'default');
            b = SY_SpreadRandomLocal(y, 100, 20, 9);
            tc.verifyNotEqual(a.stdmean, b.stdmean);
            a = NL_RQA(y, 2, 3, 10, 0.1, 2, 2, 'full', 'default');
            b = NL_RQA(y, 2, 3, 10, 0.1, 2, 2, 'full', 4);
            tc.verifyNotEqual(a.RR, b.RR); % radius from a different subsample
        end

        function fractalDimensionsReproducible(tc)
            y = Robust_finishTests.series(1200, 7);
            rng(1); a = NL_FractalDimensions(y, 2, 100, 0.2, 1, 5, {'ac', 1}, 8, {1, 3}, 'default');
            rng(2); b = NL_FractalDimensions(y, 2, 100, 0.2, 1, 5, {'ac', 1}, 8, {1, 3}, 'default');
            tc.verifyEqual(a, b);
        end

        function gaussianProcessSubsamplingReproducible(tc)
            y = Robust_finishTests.series(400, 8);
            cf = {'covSum', {'covSEiso', 'covNoise'}};
            for how = {'random_i', 'random_consec', 'random_both'}
                rng(1); a = MF_GP_Hyperparameters(y, cf, 1, 40, how{1}, 'default');
                rng(2); b = MF_GP_Hyperparameters(y, cf, 1, 40, how{1}, 'default');
                tc.verifyEqual(a, b, how{1});
            end
            rng(1); a = MF_GP_LocalPrediction(y, cf, 10, 3, 4, 'randomgap', 'default');
            rng(2); b = MF_GP_LocalPrediction(y, cf, 10, 3, 4, 'randomgap', 'default');
            tc.verifyEqual(a, b);
        end

        function delayTimeAndPreprocessingReproducible(tc)
            y = Robust_finishTests.series(800, 9);
            rng(1); a = NL_DelayTime(y, 10, 5, 'default');
            rng(2); b = NL_DelayTime(y, 10, 5, 'default');
            tc.verifyEqual(a, b);
            rng(1); a = PP_PreProcess(y, '', [], [], false, 'default');
            rng(2); b = PP_PreProcess(y, '', [], [], false, 'default');
            tc.verifyEqual(a.rmgd, b.rmgd);
            tc.verifyEqual(sort(a.rmgd), sort(BF_Random(800, 0, 'normal')), 'AbsTol', 0);
        end

        % ---------------------------------------------------------------------
        % explicit bin edges and bandwidths (no MATLAB defaults)
        % ---------------------------------------------------------------------
        function mutualInformationExplicitBins(tc)
            y = Robust_finishTests.series(2000, 11);
            a = BF_MutualInformation(y(1:end - 1), y(2:end), 'range', 'range', 10);
            tc.verifyTrue(isfinite(a) && a > 0);
            % lattice-valued data: a value exactly on an edge always goes to the same bin,
            % whatever the last bit of a rescaling does
            yl = mod(0:1999, 11)' + 0; % values 0..10, with edges of 11 bins on the lattice
            m1 = BF_MutualInformation(yl(1:end - 1), yl(2:end), 'range', 'range', 11);
            m2 = BF_MutualInformation(yl(1:end - 1) * 0.1, yl(2:end) * 0.1, 'range', 'range', 11);
            m3 = BF_MutualInformation(yl(1:end - 1) * 3, yl(2:end) * 3, 'range', 'range', 11);
            tc.verifyEqual(m1, m2, 'AbsTol', 1e-12);
            tc.verifyEqual(m1, m3, 'AbsTol', 1e-12);
            % a tiny-scale series is binned like any other (no absolute offsets)
            tc.verifyEqual(BF_MutualInformation(y(1:end - 1) * 1e-9, y(2:end) * 1e-9, 'range', 'range', 10), a, 'AbsTol', 1e-12);
            tc.verifyTrue(isnan(BF_MutualInformation(ones(100, 1), (1:100)', 'range', 'range', 10))); % constant: undefined
        end

        function distributionFitsUseExplicitDensity(tc)
            y = Robust_finishTests.series(1500, 12);
            for nb = {'sqrt', 0, 15}
                o = DN_SimpleFit(y, 'gauss1', nb{1});
                tc.verifyTrue(isfinite(o.r2) && o.r2 > 0.5, sprintf('%s', num2str(nb{1})));
                tc.verifyEqual(o, DN_SimpleFit(y, 'gauss1', nb{1}));
            end
            tc.verifyTrue(isnan(DN_SimpleFit(ones(200, 1), 'gauss1', 0))); % constant: no density
        end

        function triangularIndexAndAsymmetryUseExplicitBins(tc)
            y = Robust_finishTests.series(1000, 13);
            o = MD_RawHRVMeas(y);
            e = BF_HistEdges(y, 10);
            tc.verifyEqual(o.tri10, 1000 / max(histcounts(y, e)));
            tc.verifyEqual(MD_hrv_classic(y).tri, o.tri10);
            h = DN_HistogramAsymmetry(y, 11, false);
            tc.verifyTrue(h.modeProbPos > 0 && h.modeProbNeg > 0);
            % all values on one side of the mean: no histogram on the other side
            yo = [ones(50, 1); 2 * ones(50, 1)]; yo = (yo - mean(yo)) / std(yo);
            ho = DN_HistogramAsymmetry(yo, 11, false);
            tc.verifyEqual(ho.modeProbPos, 0.5);
        end

        function portaLevelsOnLatticeAreScaleInvariant(tc)
            [~, shuffle] = sort(BF_Random(1000, 1));
            yl = repmat((0:9)', 100, 1); yl = yl(shuffle);
            a = MD_Porta(yl, 6);
            b = MD_Porta(yl * 0.1, 6);
            tc.verifyEqual(a, b);
        end

        function walkerAndReturnTimeRun(tc)
            y = Robust_finishTests.series(1500, 14);
            w = PH_Walker(y, 'biasprop', [0.1, 0.5]);
            tc.verifyTrue(isfinite(w.sw_distdiff) && w.sw_distdiff >= 0 && w.sw_distdiff <= 2);
            r = NL_ReturnTime(y, 0.05, 100, {'ac', 1}, 500, {1, 3});
            tc.verifyTrue(isfinite(r.hhisthist) && r.maxhisthist > 0 && r.maxhisthist <= 1);
            tc.verifyEqual(r, NL_ReturnTime(y, 0.05, 100, {'ac', 1}, 500, {1, 3}));
        end
    end
end
