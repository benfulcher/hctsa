classdef Robust_degenTests < matlab.unittest.TestCase
    % Operations on degenerate input: values at the rounding level of the
    % arithmetic must come out as NaN or as a finite limiting value, never as
    % rounding noise, Inf or a complex number.

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
    end

    methods (Test)
        function dnCvNaNForZeroMean(testCase)
            rng(1);
            x = randn(500, 1);
            x = x - mean(x); % mean is zero up to rounding
            testCase.verifyTrue(isnan(DN_CV(x, 1)));
            testCase.verifyTrue(isnan(DN_CV(zscore(x), 2)));
        end

        function dnCvOrdinaryValueUnchanged(testCase)
            rng(2);
            x = 3 + randn(500, 1);
            testCase.verifyEqual(DN_CV(x, 1), std(x) / mean(x), 'AbsTol', 1e-14);
            testCase.verifyEqual(DN_CV(x, 2), (std(x) / mean(x))^2, 'AbsTol', 1e-14);
        end

        function gaussianAMIFiniteForDeterministicSeries(testCase)
            ramp = (1:300)'; % lagged copies are exactly linearly related
            ami = IN_AutoMutualInfo(ramp, 1:5, 'gaussian');
            vals = struct2array(ami);
            testCase.verifyTrue(all(isfinite(vals)) && isreal(vals));
            testCase.verifyLessThanOrEqual(max(vals), -0.5 * log(1e-12) + 1e-12);
            % an ordinary series is unaffected by the floor
            rng(3);
            y = filter(1, [1, -0.7], randn(500, 1));
            r = corr(y(1:end - 2), y(3:end));
            testCase.verifyEqual(IN_AutoMutualInfo(y, 2, 'gaussian'), -0.5 * log(1 - r^2), 'AbsTol', 0);
        end

        function spSummariesLogFieldsFiniteForPureSine(testCase)
            % a sine that fits the transform exactly has spectral bins at the rounding level
            t = (1:1024)';
            out = SP_Summaries(sin(2 * pi * t / 64), 'fft', [], [], false);
            fields = {'logstd', 'logmom3', 'logac1', 'ylogareatopeak', 'logarea_2_1', 'logiqr'};
            for i = 1:numel(fields)
                testCase.verifyTrue(isfinite(out.(fields{i})), fields{i});
            end
        end

        function varRatioTestPeriodsDoNotDependOnSaturatedPValues(testCase)
            rng(4);
            y = zscore(cumsum(randn(1000, 1)) + 5 * (1:1000)' / 1000); % strongly trending: p-values underflow
            out = SY_VarRatioTest(y, [2, 4, 6, 8, 2, 4, 6, 8], [0, 0, 0, 0, 1, 1, 1, 1]);
            [~, p, stat] = vratiotest(y, 'period', [2, 4, 6, 8, 2, 4, 6, 8], 'IID', logical([0, 0, 0, 0, 1, 1, 1, 1]));
            [~, i] = max(abs(stat));
            periods = [2, 4, 6, 8, 2, 4, 6, 8];
            testCase.verifyEqual(out.periodminpValue, periods(i));
            testCase.verifyEqual(out, SY_VarRatioTest(y, [2, 4, 6, 8, 2, 4, 6, 8], [0, 0, 0, 0, 1, 1, 1, 1]));
        end

        function nlD2NoNaNFromSinglePrecisionOverflow(testCase)
            % (needs the TISEAN binaries; a series for which the single-precision c2g overflowed)
            rng(5);
            y = zscore(cumsum(randn(1000, 1)));
            out = NL_d2(y, 1, 10, {'ac', 1});
            testCase.verifyTrue(isstruct(out));
            testCase.verifyFalse(any(isnan(struct2array(out))));
        end
    end
end
