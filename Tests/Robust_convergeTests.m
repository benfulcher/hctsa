classdef Robust_convergeTests < matlab.unittest.TestCase
    % Tests for definitions made well posed (deterministic, portable, or robust to
    % rounding) so that hctsa and pyhctsa compute the same features.

    methods (Test)
        function dvvMatchesPortableReference(tc)
            % NL_DVV draws from BF_Random: values computed by pyhctsa on the same input
            y = cumsum(BF_Random(400, 7, 'normal'));
            y = (y - mean(y)) / std(y);
            o = NL_DVV(y, 3, 50, 2, 20, 4, 0, 1);
            tc.verifyEqual(o.trend, 0.32953051863335436, 'AbsTol', 1e-10);
            tc.verifyEqual(o.max, 0.9565440085565462, 'AbsTol', 1e-10);
            tc.verifyEqual(o.meanDiff, 0.06196142307459889, 'AbsTol', 1e-10);
            tc.verifyEqual(o.rmsDiffSurr, 0.0069117290619591055, 'AbsTol', 1e-10);
            tc.verifyEqual(o.meanNormCDF, 0.3732281477819768, 'AbsTol', 1e-10);
            tc.verifyEqual(o.numZeroCrossings, 2);
        end

        function stdNthDerChangeGlobalFit(tc)
            % Differencing white noise doubles the std per order; small scale must not
            % send the exponential fit to a spurious minimum
            for scale = [1, 0.01]
                y = scale * BF_Random(1000, 11, 'normal');
                o = SY_StdNthDerChange(y, 10);
                tc.verifyGreaterThan(o.fexp_b, 0.55);
                tc.verifyLessThan(o.fexp_b, 0.75);
                tc.verifyGreaterThan(o.fexp_r2, 0.99);
            end
        end

        function compareKSFitSingularDensitiesBounded(tc)
            % Gamma/Weibull with shape < 1 and a U-shaped beta: the area between the
            % curves stays within [0, 2]
            u = BF_Random(1000, 12);
            xg = gaminv(u, 0.5, 2);
            xb = betainv(BF_Random(1000, 13), 0.4, 0.4);
            for c = {{xg, 'gamma'}, {xg, 'wbl'}, {xb, 'beta'}}
                o = DN_CompareKSFit(c{1}{1}, c{1}{2});
                tc.verifyGreaterThanOrEqual(o.adiff, 0);
                tc.verifyLessThanOrEqual(o.adiff, 2);
                tc.verifyTrue(isfinite(o.peaksepy) && isfinite(o.olapint) && isfinite(o.relent));
            end
            o = DN_CompareKSFit(BF_Random(500, 14), 'uni');
            tc.verifyTrue(isnan(o.peaksepx));
            xg(10) = 0;
            tc.verifyTrue(all(isnan(DN_CompareKSFit(xg, 'gamma'))));
        end

        function medianFilterRepeatsEnds(tc)
            % PP_Compare 'medianf3' does not depend on the offset of the data
            y = cumsum(BF_Random(500, 15, 'normal'));
            o1 = PP_Compare(y, 'medianf3');
            o2 = PP_Compare(y + 1000, 'medianf3');
            tc.verifyEqual(o1.statav2, o2.statav2, 'AbsTol', 1e-9);
            tc.verifyEqual(o1.swss2_1, o2.swss2_1, 'AbsTol', 1e-9);
        end

        function constantCurvesAreNotFitted(tc)
            % A curve of equal values (whose mean differs by rounding) gives NaN fits
            f = BF_ExpFit((1:29)', 99.887133182844 * ones(29, 1));
            tc.verifyTrue(isnan(f.r2));
            tc.verifyTrue(all(isnan(CO_AutoCorr(-0.2876717480582633 * ones(50, 1), [1, 2], 'Fourier'))));
        end

        function cepstrumIgnoresRoundingOfTheMean(tc)
            y = BF_Random(1000, 16, 'normal');
            y = (y - mean(y)) / std(y);
            o1 = SP_Cepstrum(y, 100, 4);
            o2 = SP_Cepstrum(y + 1e-15, 100, 4);
            tc.verifyEqual(o1.meanCeps, o2.meanCeps, 'AbsTol', 1e-12);
        end

        function marginalTestsLogJB(tc)
            % log of the Jarque-Bera statistic from its definition (same value in pyhctsa)
            x = [1; 2; 3; 4; 10; 2; 3; 1; 0; 5];
            tc.verifyEqual(HT_MarginalTests(x, 'logjb'), 1.5018106117098595, 'AbsTol', 1e-12);
            tc.verifyEqual(HT_MarginalTests(3*x + 7, 'logjb'), 1.5018106117098595, 'AbsTol', 1e-12);
            tc.verifyTrue(isnan(HT_MarginalTests(ones(20, 1), 'logjb')));
        end

        function saturatingFitsAreGlobal(tc)
            % SP_Summaries polysat fit: a curve that keeps rising like w^2 has b at the
            % upper limit (only a/b determined), deterministically
            y = BF_Random(1000, 17, 'normal');
            o1 = SP_Summaries(y, 'fft', [], [], false);
            o2 = SP_Summaries(y, 'fft', [], [], false);
            tc.verifyEqual(o1.fpolysat_b, o2.fpolysat_b);
            tc.verifyTrue(isfinite(o1.fpolysat_r2));
        end
    end
end
