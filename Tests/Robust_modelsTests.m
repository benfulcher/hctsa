classdef Robust_modelsTests < matlab.unittest.TestCase
    % Unit tests for model-fitting features that were made robust: bounded,
    % deterministic, and insensitive to tiny perturbations of the input.

    methods(TestClassSetup)
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

    methods(Test)

        %-------------------------------------------------------------
        % CO_PartialAutoCorr
        %-------------------------------------------------------------
        function test_CO_PartialAutoCorr_AR1(testCase)
            % AR(1): partial autocorrelation at lag 1 is the coefficient, ~0 beyond
            rng(1);
            e = randn(5000,1); y = zeros(5000,1);
            for t = 2:5000, y(t) = 0.7*y(t-1) + e(t); end
            out = CO_PartialAutoCorr(y,5);
            testCase.verifyEqual(out.pac_1, 0.7, 'AbsTol', 0.03);
            testCase.verifyLessThan(abs(out.pac_2), 0.05);
            testCase.verifyLessThan(abs(out.pac_5), 0.05);
        end

        function test_CO_PartialAutoCorr_AgreesWithOLS(testCase)
            rng(2);
            y = filter(1,[1 -0.5 0.2],randn(1000,1));
            a = CO_PartialAutoCorr(y,10); b = CO_PartialAutoCorr(y,10,'ols');
            testCase.verifyEqual(cell2mat(struct2cell(a)), cell2mat(struct2cell(b)), 'AbsTol', 0.02);
        end

        function test_CO_PartialAutoCorr_BoundedAndStableOnSinusoid(testCase)
            % An exact sinusoid is perfectly predictable by an AR(2): the partial
            % autocorrelations beyond lag 2 are ~0 and do not depend on rounding noise
            t = (1:1000)';
            y = sin(2*pi*t/50);
            rng(3);
            a = cell2mat(struct2cell(CO_PartialAutoCorr(y,20)));
            b = cell2mat(struct2cell(CO_PartialAutoCorr(y + 1e-9*randn(1000,1),20)));
            testCase.verifyLessThanOrEqual(max(abs(a)), 1);
            testCase.verifyLessThan(max(abs(a(3:end))), 0.2);
            testCase.verifyEqual(a, b, 'AbsTol', 1e-3);
        end

        function test_CO_PartialAutoCorr_Constant(testCase)
            out = CO_PartialAutoCorr(ones(100,1),3);
            testCase.verifyEqual(cell2mat(struct2cell(out)), zeros(3,1));
        end

    end
end
