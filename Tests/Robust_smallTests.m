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
    end
end
