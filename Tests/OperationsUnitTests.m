classdef OperationsUnitTests < matlab.unittest.TestCase
    % Known-answer and edge-case unit tests for individual Operations/*.m
    % functions, run in isolation (no database, no full pipeline).
    %
    % Complements BasicPipelineTests.m, which only exercises the full
    % TS_Compute pipeline end-to-end.

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
        % CO_AutoCorr
        %-------------------------------------------------------------
        function test_CO_AutoCorr_LagZeroIsOne(testCase)
            rng(1);
            y = randn(200,1);
            out = CO_AutoCorr(y,0,'Fourier');
            testCase.verifyEqual(out, 1, 'AbsTol', 1e-10, ...
                'Autocorrelation at lag 0 should be exactly 1.');
        end

        function test_CO_AutoCorr_VectorMatchesScalarAndFull(testCase)
            rng(2);
            y = randn(150,1);
            lags = [0,1,2,5,10];

            % Vector-tau call:
            outVector = CO_AutoCorr(y,lags,'Fourier');
            % Individual scalar calls:
            outScalar = arrayfun(@(tau) CO_AutoCorr(y,tau,'Fourier'), lags)';
            % Full ACF (tau = []), indexed manually:
            fullAcf = CO_AutoCorr(y,[],'Fourier');
            outFromFull = fullAcf(lags+1);

            testCase.verifyEqual(outVector, outScalar, 'AbsTol', 1e-12, ...
                'Vector-tau result should match individual scalar-tau calls.');
            testCase.verifyEqual(outVector, outFromFull, 'AbsTol', 1e-12, ...
                'Vector-tau result should match indexing into the full ACF.');
        end

        function test_CO_AutoCorr_ConstantSeriesIsNaN(testCase)
            % A constant series has zero variance, so the 0/0 normalization
            % silently yields NaN (MATLAB does not warn on array 0/0).
            y = ones(20,1);
            out = CO_AutoCorr(y,1,'Fourier');
            testCase.verifyTrue(isnan(out), 'Autocorrelation of a constant series should be NaN.');
        end

        %-------------------------------------------------------------
        % CO_FirstMin
        %-------------------------------------------------------------
        function test_CO_FirstMin_SinusoidPeriod(testCase)
            % A clean periodic signal should have its first ACF minimum
            % close to half its period.
            T = 20; % period, in samples
            N = 200;
            t = (0:N-1)';
            y = cos(2*pi*t/T);

            out = CO_FirstMin(y,'ac');

            testCase.verifyGreaterThanOrEqual(out, T/2 - 3);
            testCase.verifyLessThanOrEqual(out, T/2 + 3);
        end

        function test_CO_FirstMin_ConstantSeriesWarnsAndReturnsNaN(testCase)
            y = ones(50,1);
            lastwarn(''); % clear
            out = CO_FirstMin(y,'ac');
            [msg,~] = lastwarn();
            testCase.verifyTrue(isnan(out), 'CO_FirstMin should return NaN for a constant series.');
            testCase.verifyNotEmpty(msg, 'CO_FirstMin should warn when no minimum can be found.');
        end

        %-------------------------------------------------------------
        % CO_AutoCorrShape
        %-------------------------------------------------------------
        function test_CO_AutoCorrShape_PosDrownBasicChecks(testCase)
            rng(3);
            N = 500;
            t = (0:N-1)'/100;
            y = sin(2*pi*10*t) + 0.1*randn(N,1);

            out = CO_AutoCorrShape(y,'posDrown');

            testCase.verifyGreaterThanOrEqual(out.Nac, 1, 'Nac should be a positive lag.');
            testCase.verifyLessThanOrEqual(out.Nac, N, 'Nac cannot exceed the series length.');
            testCase.verifyFalse(isnan(out.sumacf), 'sumacf should not be NaN for a well-behaved series.');
            testCase.verifyFalse(isnan(out.meanacf), 'meanacf should not be NaN for a well-behaved series.');
        end

        function test_CO_AutoCorrShape_ConstantSeriesWarnsAndReturnsNaN(testCase)
            y = ones(50,1);
            lastwarn(''); % clear
            out = CO_AutoCorrShape(y,'posDrown');
            [msg,~] = lastwarn();
            testCase.verifyTrue(isnan(out), 'CO_AutoCorrShape should return NaN for a constant series.');
            testCase.verifyNotEmpty(msg, 'CO_AutoCorrShape should warn on an anomalous (constant) series.');
        end

        %-------------------------------------------------------------
        % SC_MMA
        %-------------------------------------------------------------
        function test_SC_MMA_RunsAndProducesFiniteOutput(testCase)
            rng(4);
            % N=3000 avoids a pre-existing edge case in SC_MMA where
            % maxScale/5 == minScale exactly (default scaleRange at N=2000),
            % which leaves an internal array empty and errors -- unrelated to
            % what's under test here.
            N = 3000;
            y = cumsum(randn(N,1));

            out = SC_MMA(y);

            testCase.verifyTrue(isfinite(out.meanHurstExponent), ...
                'meanHurstExponent should be finite for a well-behaved series.');
            testCase.verifyTrue(isfinite(out.stdHurstExponent));
            testCase.verifyGreaterThanOrEqual(out.maxHurstExponent, out.minHurstExponent);
        end

        function test_SC_MMA_ScaleBoundaryDoesNotCrash(testCase)
            % Regression test for a boundary case where the default
            % scaleRange gives maxScale/5 == minScale exactly (e.g. N=2000
            % gives scaleRange=[10,50]). This used to make an internal
            % zero-step colon range (minScale:0:(maxScale/5)) silently
            % return empty, crashing later at find(hqs==max(hqs(:)),1).
            rng(5);
            N = 2000;
            y = cumsum(randn(N,1));

            out = SC_MMA(y);

            testCase.verifyTrue(isfinite(out.meanHurstExponent));
        end

        %-------------------------------------------------------------
        % DN_Mean
        %-------------------------------------------------------------
        function test_DN_Mean_KnownAnswers(testCase)
            y = (1:5)';

            testCase.verifyEqual(DN_Mean(y,'arithmetic'), 3, 'AbsTol', 1e-10);
            testCase.verifyEqual(DN_Mean(y,'median'), 3, 'AbsTol', 1e-10);
            testCase.verifyEqual(DN_Mean(y,'geom'), geomean(y), 'AbsTol', 1e-10);
            testCase.verifyEqual(DN_Mean(y,'harm'), harmmean(y), 'AbsTol', 1e-10);
            testCase.verifyEqual(DN_Mean(y,'rms'), sqrt(sum(y.^2)/length(y)), 'AbsTol', 1e-10);
        end

        function test_DN_Mean_UnknownTypeErrors(testCase)
            y = (1:5)';
            threw = false;
            try
                DN_Mean(y,'not_a_real_type');
            catch
                threw = true;
            end
            testCase.verifyTrue(threw, 'DN_Mean should error for an unrecognized meanType.');
        end

    end
end
