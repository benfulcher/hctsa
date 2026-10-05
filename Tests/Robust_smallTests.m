classdef Robust_smallTests < matlab.unittest.TestCase
    % Tests of MF_CompareAR's bounded loss ratios and of the rounding-level guards in CO_FallingSticks.

    methods (Test)
        function compareARRatiosAreBounded(tc)
            rng(1);
            t = (1:600)';
            ys = {randn(600, 1), cumsum(randn(600, 1)), sin(0.1*t) + 0.5*randn(600, 1), sin(0.1*t)};
            for i = 1:numel(ys)
                for testHow = {0.5, 'all'}
                    o = MF_CompareAR(zscore(ys{i}), 1:10, testHow{1});
                    tc.verifyGreaterThanOrEqual(o.propgain1min, 0);
                    tc.verifyLessThanOrEqual(o.propgain1min, 1);
                    tc.verifyGreaterThan(o.medonmax, 0);
                    tc.verifyLessThanOrEqual(o.medonmax, 1);
                    tc.verifyFalse(isfield(o, 'firstonmin') || isfield(o, 'maxonmed'));
                end
            end
        end

        function compareARExactFitGivesLimit(tc)
            % a noiseless two-sine series is predicted exactly by a higher order: AR(1) error is all removed
            t = (1:800)';
            o = MF_CompareAR(zscore(sin(0.1*t) + sin(0.37*t)), 1:10, 'all');
            tc.verifyGreaterThan(o.propgain1min, 1 - 1e-9);
            tc.verifyLessThan(o.medonmax, 1e-9);
        end

        function compareARWhiteNoiseGainsLittle(tc)
            rng(2);
            o = MF_CompareAR(randn(1000, 1), 1:10, 0.5);
            tc.verifyLessThan(o.propgain1min, 0.1);
            tc.verifyGreaterThan(o.medonmax, 0.9);
        end

        function fallingSticksFlatBranchPersistenceIsNaN(tc)
            % all negative sticks are shorter than their spacing and fall flat (pi/2): the angles are
            % constant, but their mean rounds, so std(angles) is not exactly zero
            rng(3);
            y = [0.3*rand(300, 1); -0.3*rand(700, 1)];
            y = y(randperm(1000));
            o = CO_FallingSticks(y);
            tc.verifyEqual(o.propFlat_n, 1);
            tc.verifyTrue(isnan(o.tau_n) && isnan(o.ac1_n));
            tc.verifyTrue(isnan(o.tau_p) && isnan(o.ac1_p)); % so is the positive branch
            tc.verifyEqual(o.std_all, 0);
            tc.verifyTrue(isnan(o.skewness_all) && isnan(o.kurtosis_all));
            tc.verifyEqual(o.q10_all, pi/2);
        end

        function fallingSticksOrdinarySeriesUnaffected(tc)
            rng(4);
            o = CO_FallingSticks(zscore(randn(1000, 1)));
            tc.verifyFalse(any(isnan([o.tau_p o.tau_n o.ac1_p o.ac1_n o.skewness_all o.kurtosis_all])));
            tc.verifyGreaterThan(o.std_all, 0.05);
        end
    end
end
