classdef Robust_binsTests < matlab.unittest.TestCase
    % Tests of the explicit-edge histogram and explicit-bandwidth kernel density
    % helpers and of the operations built on them: bin assignment must not depend
    % on the last bit of the data or of the edges, and the definitions must be
    % bounded and deterministic.

    methods (Test)
        function histEdgesLatticeIsStable(tc)
            % Lattice-valued data sit exactly on ideal edges: a one-ulp change must not move them.
            y = repmat((0:10)'/10, 30, 1);
            e = BF_HistEdges(y, 10);
            c1 = histcounts(y, e);
            c2 = histcounts(y*(1 + 4*eps), e);
            c3 = histcounts(y*(1 - 4*eps), e);
            tc.verifyEqual(c1, c2);
            tc.verifyEqual(c1, c3);
            tc.verifyEqual(sum(c1), numel(y));
        end

        function histEdgesRules(tc)
            rng(1); y = randn(1000, 1);
            tc.verifyEqual(numel(BF_HistEdges(y, 'sqrt')) - 1, ceil(sqrt(1000)));
            tc.verifyEqual(numel(BF_HistEdges(y, 'sturges')) - 1, ceil(log2(1000) + 1));
            tc.verifyGreaterThanOrEqual(numel(BF_HistEdges(y, 'auto')), numel(BF_HistEdges(y, 'sturges')));
            e = BF_HistEdges(y, 7, [-3 3]);
            tc.verifyEqual(numel(e), 8);
            tc.verifyEqual(e(2) - e(1), 6/7, 'AbsTol', 1e-5);
        end

        function quantileEdgesAreEquiprobable(tc)
            rng(2); y = randn(1000, 1);
            c = histcounts(y, BF_QuantileEdges(y, 10));
            tc.verifyEqual(c, 100*ones(1, 10), 'AbsTol', 1);
        end

        function ksDensityMatchesFormula(tc)
            rng(3); y = randn(200, 1);
            [f, xi, h] = BF_KSDensity(y);
            tc.verifyEqual(h, median(abs(y - median(y)))/0.6745 * (4/(3*200))^(1/5), 'RelTol', 1e-12);
            tc.verifyEqual(sum(f)*(xi(2) - xi(1)), 1, 'AbsTol', 1e-3);
            tc.verifyTrue(all(isnan(BF_KSDensity(ones(10, 1), 0))));
        end

        function halfSampleMode(tc)
            tc.verifyEqual(BF_HalfSampleMode([1 2 2 2 3 10 11]'), 2);
            tc.verifyEqual(BF_HalfSampleMode(5), 5);
            rng(4); y = randn(2000, 1) + 3;
            tc.verifyEqual(BF_HalfSampleMode(y), 3, 'AbsTol', 0.3);
            tc.verifyEqual(DN_CustomSkewness(y, 'pearsonMode'), DN_CustomSkewness(y + 1e-9*randn(size(y)), 'pearsonMode'), 'AbsTol', 1e-3);
        end

        function stickAnglesAndLocalDistributionsBounded(tc)
            rng(5); y = zscore(cumsum(randn(1000, 1)) + randn(1000, 1));
            o = CO_StickAngles(y);
            tc.verifyGreaterThanOrEqual(o.pnsumabsdiff, 0); tc.verifyLessThanOrEqual(o.pnsumabsdiff, 2);
            tc.verifyGreaterThanOrEqual(o.symks_p, 0); tc.verifyLessThanOrEqual(o.symks_p, 1);
            tc.verifyGreaterThanOrEqual(o.symks_n, 0); tc.verifyLessThanOrEqual(o.symks_n, 1);
            d = SY_LocalDistributions(y, 5, 'par');
            tc.verifyGreaterThanOrEqual(d.meandiv, 0); tc.verifyLessThanOrEqual(d.meandiv, 1);
            % quantized data (many tied values) still give a defined answer
            yq = round(2*y);
            tc.verifyTrue(isfinite(SY_LocalDistributions(yq, 5, 'par').meandiv));
            tc.verifyTrue(isfinite(CO_StickAngles(zscore(yq)).pnsumabsdiff));
        end

        function distributionEntropyRulesAndDegenerate(tc)
            rng(6); y = randn(1000, 1);
            tc.verifyEqual(EN_DistributionEntropy(y, 'hist', 'sqrt'), EN_DistributionEntropy(y, 'hist', ceil(sqrt(1000))), 'AbsTol', 1e-12);
            tc.verifyTrue(isnan(EN_DistributionEntropy(ones(100, 1), 'hist', 10)));
        end

        function amiEdgesAndPeakCount(tc)
            % tied, lattice-valued data: bin edges at quantiles sit on data values
            rng(7); y = zscore(round(3*randn(600, 1)));
            a = CO_HistogramAMI(y, 1, 'quantiles', 5);
            b = CO_HistogramAMI(y*(1 + 4*eps), 1, 'quantiles', 5);
            tc.verifyEqual(a, b, 'AbsTol', 1e-12);
            o = CO_CompareMinAMI(zscore(cumsum(randn(500, 1))), 'even', 2:30);
            tc.verifyGreaterThanOrEqual(o.nprompeaks, 0);
        end

        function embed2HistogramsOnEdges(tc)
            % tied values give angles of exactly 0 and +/-pi/2, which sit on histogram edges
            y = zscore(repmat([0 1 1 2 0 0 2 1]', 100, 1));
            o = CO_Embed2(y, 1);
            tc.verifyEqual(o.histent, CO_Embed2(y*(1 + 4*eps), 1).histent, 'AbsTol', 1e-9);
            s = CO_Embed2_Shapes(y, 1, 'circle', 0.5, 0);
            tc.verifyTrue(isfinite(s.hist_ent) && isfinite(s.mode));
        end

        function slidingWindowSpectralEntropy(tc)
            rng(8); y = zscore(randn(1000, 1));
            v = SY_SlidingWindow(y, 'specen', 'std', 5, 1);
            tc.verifyTrue(isfinite(v));
            tc.verifyEqual(v, SY_SlidingWindow(y, 'specen', 'std', 5, 1));
        end

        function kernelSmoothDegenerate(tc)
            tc.verifyTrue(isnan(DN_FitKernelSmooth(ones(50, 1))));
            rng(9); o = DN_FitKernelSmooth(randn(500, 1));
            tc.verifyTrue(isfinite(o.entropy) && isfinite(o.asym));
        end
    end
end
