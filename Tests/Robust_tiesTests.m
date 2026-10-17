classdef Robust_tiesTests < matlab.unittest.TestCase
    % Exact repeats and tied distances in the nonlinear operations that search
    % for neighbors. Distances within 1e-12 of the range of the data are equal,
    % a pair closer than that is a repeated point, and equally near neighbors
    % are treated alike. So on series with exact repeats (a periodic orbit
    % reached after a transient, an integer-valued series) the outputs are
    % identical, up to rounding, for y, 3.7*y - 12 and -y (which z-score to the
    % same series up to rounding and sign); quantities that cannot be estimated
    % from such series are NaN rather than numbers set by rounding; and on
    % continuous data the outputs are unchanged. The series are the same as in
    % pyhctsa's tests/test_nl_ties.py.

    properties (TestParameter)
        opName = {'FractalDimensions', 'LargestLyap', 'ReturnTime', 'DVV', 'EVTLocalDim', 'LyapSpec', 'c1'};
        tieSeries = {'orbit', 'integer', 'sin7sin20'};
    end

    methods (TestClassSetup)
        function setPaths(~)
            here = fileparts(fileparts(mfilename('fullpath')));
            setenv('PATH', [fullfile(here, 'Toolboxes', 'Tisean_3.0.1', 'bin'), ':', getenv('PATH')]);
        end
    end

    methods (Test)
        function tiesDoNotDependOnRounding(tc, opName, tieSeries)
            % (NL_LyapSpec adds a dither of fixed sign pattern, so only the affine
            % map; the fits of NL_FractalDimensions are found to the optimizer's
            % tolerance, and lyap_r writes its curve to 6 significant digits, so
            % those two agree to 1e-4)
            x = testSeries(tieSeries);
            ref = runOp(opName, x);
            relTol = 1e-6;
            if ismember(opName, {'FractalDimensions', 'LargestLyap'}), relTol = 1e-4; end
            verifySame(tc, runOp(opName, 3.7*x - 12), ref, relTol);
            if ~strcmp(opName, 'LyapSpec')
                verifySame(tc, runOp(opName, -x), ref, relTol);
            end
        end

        function fractalDimensionsNotIdentifiedOnPeriodicOrbit(tc)
            % the k-th neighbor distances are all equal or all zero: the moments
            % do not grow with k (the error has its minimum at the upper end of
            % the search interval) or are zero (the error is flat), and D is not
            % identified
            for s = {'sin7', 'sin7sin20'}
                out = runOp('FractalDimensions', testSeries(s{1}));
                tc.verifyTrue(all(isnan(fieldValues(out))), s{1});
            end
        end

        function returnTimeShortLagSharesAreBounded(tc)
            % the profile of this exactly periodic series is zero over the second
            % half of the lags: all of its mass and variation are at short lags,
            % a share of 1 (never infinite)
            out = runOp('ReturnTime', testSeries('sin7sin20'));
            tc.verifyEqual([out.fracmeanrt, out.fracsdrt], [1, 1]);
            tc.verifyFalse(any(isinf(fieldValues(out))));
        end

        function dvvTargetVarianceIsNeverNegative(tc)
            for s = {'sin7', 'orbit', 'integer', 'sin7sin20'}
                out = runOp('DVV', testSeries(s{1}));
                tc.verifyGreaterThanOrEqual(out.min, 0, s{1});
            end
            out = runOp('DVV', testSeries('sin7'));
            tc.verifyEqual(out.min, 0); % targets of equal delay vectors are equal
        end

        function evtLocalDimIsNaNForDiscreteDistances(tc)
            % the distances to a pole of an integer-valued series take a few
            % values: no continuous tail
            out = runOp('EVTLocalDim', testSeries('integer'));
            tc.verifyTrue(isnan(out.meanLocalDim));
            tc.verifyLessThan(out.propValidPoles, 0.5);
        end

        function c1IsNaNOnPeriodicOrbits(tc)
            % at every mass some reference point coincides with its k-th neighbor
            for s = {'sin7', 'orbit'}
                out = runOp('c1', testSeries(s{1}));
                tc.verifyTrue(all(isnan(fieldValues(out))), s{1});
            end
        end

        function largestLyapOnPeriodicOrbitIsBounded(tc)
            % repeated points are not neighbors: the nearest neighbor is another
            % state of the orbit, and the divergence curve is bounded rather than
            % built from coincident pairs
            out = runOp('LargestLyap', testSeries('sin7sin20'));
            tc.verifyTrue(isfinite(out.maxp));
            tc.verifyLessThan(abs(out.maxp), 1);
        end

        function knownValuesOnQuantizedChaoticMaps(tc)
            % recording a chaotic map to two decimals creates many ties and
            % repeated points, but must not spoil the estimates: the largest
            % exponent of the logistic map at r = 4 is ln 2 and that of the Henon
            % map about 0.42, and the information dimension of the Henon
            % attractor about 1.2
            x = testSeries('henonLong');
            xq = round(x, 2);
            out = runOp('LargestLyap', logisticMap(4, 0.3141, 200, 5000));
            tc.verifyEqual(out.vse_gradient, log(2), 'AbsTol', 0.03);
            out = runOp('LargestLyap', xq);
            tc.verifyEqual(out.vse_gradient, 0.42, 'AbsTol', 0.05);
            out = runOp('c1', xq);
            tc.verifyEqual(out.mediand, 1.22, 'AbsTol', 0.1);
            out = runOp('LyapSpec', xq);
            tc.verifyEqual(out.LE1, 0.42, 'AbsTol', 0.05);
            a = runOp('FractalDimensions', x); b = runOp('FractalDimensions', xq);
            tc.verifyEqual(b.meanDq, a.meanDq, 'AbsTol', 0.1);
        end

        function continuousSeriesUnchanged(tc)
            % values on continuous series, which have no repeated points or ties
            E = controlValues();
            for i = 1:size(E, 1)
                out = runOp(E{i, 2}, testSeries(E{i, 1}));
                f = fieldnames(E{i, 3});
                for j = 1:numel(f)
                    tc.verifyEqual(out.(f{j}), E{i, 3}.(f{j}), 'RelTol', 1e-9, 'AbsTol', 1e-12, ...
                        sprintf('%s, %s.%s', E{i, 1}, E{i, 2}, f{j}));
                end
            end
        end
    end
end

% ------------------------------------------------------------------------------
function out = runOp(name, x)
    x = zscore(x);
    switch name
        case 'FractalDimensions'
            out = NL_FractalDimensions(x, 2, 10, -1, 1, 5, {'ac', 1}, 32, {1, 5}, 'default');
        case 'LargestLyap'
            out = NL_LargestLyap(x, -1, {'ac1e', 30}, {'ac1e', 1}, 3, {1, 4});
        case 'ReturnTime'
            out = NL_ReturnTime(x, 0.02, 100, {'ac', 1}, 500, {1, 8});
        case 'DVV'
            out = NL_DVV(x, 3, 100, 2, 50, 10, 'default');
        case 'EVTLocalDim'
            out = NL_EVTLocalDim(x, 'ac1e', 3, 0.98, {'ac', 1}, 200, 5, 'full', 'default');
        case 'LyapSpec'
            out = NL_LyapSpec(x, 1, 3, 30, 'full', {'ac', 1});
        case 'c1'
            out = NL_c1(x, 1, [2, 6], 25, 0.1);
    end
end

function verifySame(tc, a, b, relTol)
    % the same NaN pattern, and values equal up to rounding
    if ~isstruct(a), a = struct('value', a); end
    if ~isstruct(b), b = struct('value', b); end
    va = fieldValues(a); vb = fieldValues(b);
    if isscalar(va) && isnan(va), va = nan(size(vb)); end
    if isscalar(vb) && isnan(vb), vb = nan(size(va)); end
    tc.verifyEqual(isnan(va), isnan(vb));
    ok = ~isnan(va) & ~isnan(vb);
    tc.verifyEqual(va(ok), vb(ok), 'RelTol', relTol, 'AbsTol', 1e-9);
end

function x = testSeries(name)
    switch name
        case 'orbit' % a period-3 orbit after a short transient
            x = logisticMap(3.83, 0.3, 0, 1200);
        case 'integer' % 9 integer levels
            x = floor(8 * logisticMap(3.9, 0.4, 100, 1500));
        case 'sin7' % exactly periodic
            x = repmat(sin(2*pi*(0:6)'/7 + 0.3), 150, 1);
        case 'sin7sin20' % exactly periodic, period 140
            t7 = sin(2*pi*(0:6)'/7 + 0.3); t20 = 0.5*sin(2*pi*(0:19)'/20 + 0.9);
            i = (0:1399)';
            x = t7(mod(i, 7) + 1) + t20(mod(i, 20) + 1);
        case {'henon', 'henonLong'}
            if strcmp(name, 'henon'), n = 1500; else, n = 5000; end
            a = 1.4; b = 0.3; u = 0.1; v = 0.1;
            for i = 1:100, w = 1 - a*u*u + v; v = b*u; u = w; end
            x = zeros(n, 1);
            for i = 1:n, x(i) = u; w = 1 - a*u*u + v; v = b*u; u = w; end
        case 'lorenz' % sigma = 10, rho = 28, beta = 8/3; RK4, step 0.01, every 2 steps
            s = [1; 1; 1]; h = 0.01;
            f = @(s) [10*(s(2)-s(1)); s(1)*(28-s(3))-s(2); s(1)*s(2)-8/3*s(3)];
            for i = 1:5000, s = rk4(f, s, h); end
            x = zeros(1500, 1);
            for i = 1:1500
                s = rk4(f, rk4(f, s, h), h);
                x(i) = s(1);
            end
        case 'ar'
            e = BF_Random(1600, 7, 'normal');
            x = zeros(1600, 1);
            for i = 2:1600, x(i) = 0.8*x(i-1) + e(i); end
            x = x(101:end);
    end
end

function v = fieldValues(s)
    % the values of an output structure (a scalar output, such as NaN, as it is)
    if isstruct(s), v = cell2mat(struct2cell(s))'; else, v = s; end
end

function x = logisticMap(r, x0, burn, n)
    x = zeros(n, 1); v = x0;
    for i = 1:burn, v = r*v*(1 - v); end
    for i = 1:n, x(i) = v; v = r*v*(1 - v); end
end

function s = rk4(f, s, h)
    k1 = f(s); k2 = f(s + h/2*k1); k3 = f(s + h/2*k2); k4 = f(s + h*k3);
    s = s + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

function E = controlValues()
    E = {
        'henon', 'FractalDimensions', struct('maxDq', 1.5016375196984, 'meanDq', 1.38521900424072, 'linfit_a', -0.0758464254456222);
        'henon', 'LargestLyap', struct('p2', 0.233489000000001, 'maxp', 5.1667239, 'vse_gradient', 0.422918766666667);
        'henon', 'ReturnTime', struct('std', 0.471509540881739, 'fracmeanrt', 0.4590929132545169, 'hhist', 4.4894830130248);
        'henon', 'DVV', struct('min', 0.00745023096111742, 'mean', 0.712964878798353, 'meanNormCDF', 0.0712135111771347);
        'henon', 'EVTLocalDim', struct('meanLocalDim', 1.41606998272907, 'meanTheta', 0.977744966442953);
        'henon', 'LyapSpec', struct('LE1', 0.4273765, 'LE2', -0.2992378, 'KYdim', 2.07522870015182);
        'henon', 'c1', struct('bestestd', 1.22310799, 'mediand', 1.24000078909091);
        'lorenz', 'FractalDimensions', struct('maxDq', 1.9751703682233, 'meanDq', 1.89145847456705, 'linfit_a', 0.082605593889623);
        'lorenz', 'LargestLyap', struct('p2', 0.105110000000001, 'maxp', 3.945116, 'vse_gradient', 0.0149060371019573);
        'lorenz', 'ReturnTime', struct('std', 1.14002587757889, 'fracmeanrt', 0.49613360516532407, 'hhist', 3.7700599732039);
        'lorenz', 'DVV', struct('min', 0.000728125801039375, 'mean', 0.451756458674478, 'meanNormCDF', 0.170153772311658);
        'lorenz', 'EVTLocalDim', struct('meanLocalDim', 2.20771651829273, 'meanTheta', 0.302932913651838);
        'lorenz', 'LyapSpec', struct('LE1', 0.01803355, 'LE2', -0.05220342, 'KYdim', 1.34544767373479);
        'lorenz', 'c1', struct('bestestd', 2.062404215, 'mediand', 2.00544669333333);
        'ar', 'FractalDimensions', struct('maxDq', 5.03606326170433, 'meanDq', 4.97704712814451, 'linfit_a', 0.110370161842856);
        'ar', 'LargestLyap', struct('p2', 0.875545, 'maxp', 2.165377, 'vse_gradient', 0.268768);
        'ar', 'ReturnTime', struct('std', 0.346663014385052, 'fracmeanrt', 0.5353737761818226, 'hhist', 4.55342960769032);
        'ar', 'DVV', struct('min', 0.371304904635249, 'mean', 0.697908079471809, 'meanNormCDF', 0.817785250649406);
        'ar', 'EVTLocalDim', struct('meanLocalDim', 3.32596973866923, 'meanTheta', 0.737900979387082);
        'ar', 'LyapSpec', struct('LE1', 0.01153678, 'LE2', -0.3547304, 'KYdim', 1.03252267073812);
        'ar', 'c1', struct('bestestd', 2.054448389, 'mediand', 3.94965887166667);
    };
end
