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

        function clearSharedCaches(testCase)
            % CO_AutoCorr.m and CO_FirstMin.m cache their results across
            % calls (persistent variables) to avoid redundant recomputation
            % across operations. Clear them here so tests are hermetic and
            % order-independent regardless of how many times this test class
            % runs within the same live MATLAB session (a cache hit skips not
            % just the recomputation but also any warning() side effect the
            % computation would have raised, which a couple of tests below
            % check for directly).
            clear CO_AutoCorr CO_FirstMin NL_crptool_fnn BF_CheckToolbox CO_glscf
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
        % CO_AutoCorr caching
        %-------------------------------------------------------------
        function test_CO_AutoCorr_CacheDoesNotCrossContaminate(testCase)
            % CO_AutoCorr caches the last few computed ACFs (keyed on an
            % exact, NaN-safe content match of y) so repeated calls for the
            % same series can skip recomputing the FFT. Confirm interleaved
            % calls with different series don't get confused with each
            % other's cached values, against an independent reference
            % computation (the original, uncached algorithm).
            rng(13);
            y1 = randn(200,1);
            y2 = randn(250,1);

            out1a = CO_AutoCorr(y1,[],'Fourier');
            out2 = CO_AutoCorr(y2,[],'Fourier');
            out1b = CO_AutoCorr(y1,[],'Fourier'); % should hit the cache

            testCase.verifyEqual(out1a, referenceFourierACF(y1), 'AbsTol', 1e-10);
            testCase.verifyEqual(out2, referenceFourierACF(y2), 'AbsTol', 1e-10);
            testCase.verifyEqual(out1b, referenceFourierACF(y1), 'AbsTol', 1e-10);
            testCase.verifyEqual(out1a, out1b, 'AbsTol', 0, ...
                'Repeated calls for the same series should give identical results.');
        end

        function test_CO_AutoCorr_CacheHandlesConstantSeries(testCase)
            % A constant series produces all-NaN ACF entries (0/0
            % normalization). isequaln (not isequal) is used for the
            % cache-match check specifically so a repeated call with the same
            % NaN-containing series still hits the cache rather than being
            % treated as perpetually distinct.
            y = ones(30,1);
            out1 = CO_AutoCorr(y,[],'Fourier');
            out2 = CO_AutoCorr(y,[],'Fourier'); % should hit the cache
            testCase.verifyTrue(all(isnan(out1)));
            testCase.verifyEqual(out1, out2); % verifyEqual treats NaN==NaN as equal
        end

        function test_CO_AutoCorr_CacheEvictionStillCorrect(testCase)
            % Push more distinct series through than the cache holds, then
            % confirm a call for the very first (now-evicted) series still
            % recomputes correctly rather than returning stale/wrong data.
            rng(14);
            N = 100;
            firstSeries = randn(N,1);
            out1a = CO_AutoCorr(firstSeries,[],'Fourier');
            for k = 1:10
                CO_AutoCorr(randn(N,1),[],'Fourier'); % evict firstSeries from the cache
            end
            out1b = CO_AutoCorr(firstSeries,[],'Fourier'); % must recompute correctly

            testCase.verifyEqual(out1a, out1b, 'AbsTol', 1e-10);
            testCase.verifyEqual(out1a, referenceFourierACF(firstSeries), 'AbsTol', 1e-10);
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

        function test_CO_FirstMin_CacheDoesNotCrossContaminate(testCase)
            % CO_FirstMin caches its resolved output (keyed on an exact,
            % NaN-safe match of y, minWhat, extraParam, minNotMax). Confirm
            % interleaved calls with different series don't get confused with
            % each other's cached values, against an independent reference
            % reimplementation of the 'ac'-branch algorithm (no caching).
            rng(15);
            y1 = randn(150,1);
            y2 = randn(180,1);

            out1a = CO_FirstMin(y1,'ac');
            out2 = CO_FirstMin(y2,'ac');
            out1b = CO_FirstMin(y1,'ac'); % should hit the cache

            testCase.verifyEqual(out1a, referenceFirstMinAC(y1));
            testCase.verifyEqual(out2, referenceFirstMinAC(y2));
            testCase.verifyEqual(out1b, referenceFirstMinAC(y1), ...
                'A cache hit for y1 must not be contaminated by the intervening call for y2.');
        end

        function test_CO_FirstMin_CacheEvictionStillCorrect(testCase)
            % Push more distinct (series,method) combinations through than
            % the cache holds, then confirm a call for the very first
            % (now-evicted) combination still recomputes correctly.
            rng(16);
            N = 150;
            firstSeries = randn(N,1);
            out1a = CO_FirstMin(firstSeries,'ac');
            for k = 1:10
                CO_FirstMin(randn(N,1),'ac'); % evict firstSeries from the cache
            end
            out1b = CO_FirstMin(firstSeries,'ac'); % must recompute correctly

            testCase.verifyEqual(out1a, out1b);
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
        % SY_RangeEvolve
        %-------------------------------------------------------------
        function test_SY_RangeEvolve_MatchesBruteForcePrefixRange(testCase)
            % SY_RangeEvolve used to compute cums(i) = range(y(1:i)) via a
            % per-i call to range() on the growing prefix (O(N^2)); it now
            % uses cummax(y)-cummin(y) (O(N)). Cross-check derived output
            % fields against an independent brute-force reimplementation of
            % cums to confirm the values are unchanged.
            rng(7);
            N = 1500;
            y = randn(N,1);

            cumsExpected = zeros(N,1);
            for i = 1:N
                cumsExpected(i) = range(y(1:i));
            end
            fullrExpected = range(y);

            out = SY_RangeEvolve(y);

            testCase.verifyEqual(out.totnuq, length(unique(cumsExpected)), 'AbsTol', 0);
            testCase.verifyEqual(out.p50, cumsExpected(ceil(N*0.5))/fullrExpected, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.l100, cumsExpected(100)/fullrExpected, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.l1000, cumsExpected(1000)/fullrExpected, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.nuqp10, length(unique(cumsExpected(1:floor(N*0.1))))/out.totnuq, 'AbsTol', 1e-12);
        end

        %-------------------------------------------------------------
        % WL_DetailCoeffs
        %-------------------------------------------------------------
        function test_WL_DetailCoeffs_MatchesPerLevelDecomposition(testCase)
            % WL_DetailCoeffs used to call wavedec(y,level,wname) fresh at
            % every level (redoing the full decomposition from scratch each
            % time); it now decomposes once at maxlevel and extracts each
            % level's detail via wrcoef from that single decomposition.
            % Cross-check against an independent per-level reimplementation
            % (the original approach) to confirm identical output.
            rng(8);
            N = 2000;
            y = randn(N,1);
            wname = 'db3';
            maxlevel = 5;

            means = zeros(maxlevel,1);
            medians = zeros(maxlevel,1);
            maxs = zeros(maxlevel,1);
            for k = 1:maxlevel
                [c,l] = wavedec(y,k,wname);
                det = wrcoef('d',c,l,wname,k);
                means(k) = mean(abs(det));
                medians(k) = median(abs(det));
                maxs(k) = max(abs(det));
            end
            means_s = sort(means,'descend');
            medians_s = sort(medians,'descend');
            maxs_s = sort(maxs,'descend');

            out = WL_DetailCoeffs(y,wname,maxlevel);

            testCase.verifyEqual(out.max_mean, means_s(1), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.max_median, medians_s(1), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.max_max, maxs_s(1), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.std_mean, std(means), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.std_median, std(medians), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.std_max, std(maxs), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.max1on2_mean, means_s(1)/means_s(2), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.max1on2_median, medians_s(1)/medians_s(2), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.max1on2_max, maxs_s(1)/maxs_s(2), 'AbsTol', 1e-10);
            r = corrcoef(maxs,medians);
            testCase.verifyEqual(out.corrcoef_max_medians, r(1,2), 'AbsTol', 1e-10);
        end

        %-------------------------------------------------------------
        % NW_VisibilityGraph
        %-------------------------------------------------------------
        function test_NW_VisibilityGraph_MatchesDenseBruteForce(testCase)
            % NW_VisibilityGraph's 'horiz' method used to build an N x N dense
            % adjacency matrix (mostly zeros, since the graph has only O(N)
            % edges); it now builds a sparse matrix from accumulated edge
            % indices. Cross-check the resulting degree distribution (and a
            % couple of downstream stats) against an independent dense
            % reimplementation of the original algorithm.
            rng(9);
            N = 300;
            y = randn(N,1);
            y = y - min(y);

            % Brute-force dense reimplementation (the original algorithm):
            Aexpected = zeros(N);
            yr = flipud(y);
            for i = 1:N
                if i < N
                    nAhead = find(y(i+1:end) > y(i),1,'first');
                    if ~isempty(nAhead)
                        Aexpected(i,i+nAhead) = 1;
                    end
                end
                if i > 1
                    nBack = find(yr(N-i+2:end) > yr(N-i+1),1,'first');
                    if ~isempty(nBack)
                        Aexpected(i-nBack,i) = 1;
                    end
                end
            end
            At = Aexpected';
            lowerT = logical(tril(ones(size(Aexpected))));
            Aexpected(lowerT) = At(lowerT);
            kExpected = full(sum(Aexpected));

            out = NW_VisibilityGraph(y,'horiz');

            testCase.verifyEqual(out.meank, mean(kExpected), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.maxk, max(kExpected), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.mink, min(kExpected), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.stdk, std(kExpected), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.modek, mode(kExpected), 'AbsTol', 1e-10);
        end

        %-------------------------------------------------------------
        % CO_Embed2_Shapes
        %-------------------------------------------------------------
        function test_CO_Embed2_Shapes_MatchesExplicitOnesBroadcast(testCase)
            % CO_Embed2_Shapes used m - ones(N,1)*m(i,:) to broadcast-subtract
            % a row from every row of m; it now relies on implicit
            % broadcasting (m - m(i,:)), which should give an identical
            % result. Cross-check the 'circle' counts against the original
            % explicit-ones approach.
            rng(10);
            N = 400;
            y = randn(N,1);
            tau = 3;
            r = 1;

            m = [y(1:end-tau), y(1+tau:end)];
            Nm = length(m);
            countsExpected = zeros(Nm,1);
            for i = 1:Nm
                m_c = m - ones(Nm,1)*m(i,:);
                m_c_d = sum(m_c.^2,2);
                countsExpected(i) = sum(m_c_d <= r^2);
            end
            countsExpected = countsExpected - 1;

            out = CO_Embed2_Shapes(y,tau,'circle',r);

            testCase.verifyEqual(out.mean, mean(countsExpected), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.max, max(countsExpected), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.median, median(countsExpected), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.std, std(countsExpected), 'AbsTol', 1e-10);
        end

        %-------------------------------------------------------------
        % EN_ApEn
        %-------------------------------------------------------------
        function test_EN_ApEn_MatchesBruteForceLoop(testCase)
            % EN_ApEn used explicit for-loops to (a) build the Hankel-style
            % embedding matrix x and (b) broadcast each row x(i,:) across an
            % ax matrix before computing distances; both are now vectorized
            % (implicit broadcasting / direct indexing). Cross-check against
            % the original brute-force implementation.
            rng(11);
            N = 300;
            y = randn(N,1);
            mnom = 2;
            rth = 0.2;

            r = rth*std(y);
            phiExpected = zeros(2,1);
            for k = 1:2
                m = mnom+k-1;
                C = zeros(N-m+1,1);
                x = zeros(N-m+1,m);
                for i = 1:N-m+1
                    x(i,:) = y(i:i+m-1);
                end
                ax = ones(N-m+1,m);
                for i = 1:N-m+1
                    for j = 1:m
                        ax(:,j) = x(i,j);
                    end
                    d = abs(x-ax);
                    if m > 1
                        d = max(d,[],2)';
                    end
                    dr = (d<=r);
                    C(i) = sum(dr)/(N-m+1);
                end
                phiExpected(k) = mean(log(C));
            end
            expected = phiExpected(1) - phiExpected(2);

            actual = EN_ApEn(y,mnom,rth);

            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-10);
        end

        %-------------------------------------------------------------
        % EN_PermEn
        %-------------------------------------------------------------
        function test_EN_PermEn_MatchesBruteForceSearch(testCase)
            % EN_PermEn used to linearly scan permList (all m! permutations)
            % for each embedding vector to find its matching permutation; it
            % now uses a hash lookup built once from permList's own rows.
            % Cross-check against the original brute-force linear search.
            rng(12);
            N = 500;
            y = randn(N,1);
            m = 4;
            tau = 1;

            x = BF_Embed(y,tau,m,0);
            Nx = size(x,1);
            permList = perms(1:m);
            numPerms = length(permList);
            countPermsExpected = zeros(numPerms,1);
            for j = 1:Nx
                [~,ix] = sort(x(j,:));
                for k = 1:numPerms
                    if all(permList(k,:)-ix == 0)
                        countPermsExpected(k) = countPermsExpected(k) + 1;
                        break
                    end
                end
            end
            pExpected = countPermsExpected/Nx;
            p0Expected = pExpected(pExpected>0);
            permEnExpected = -sum(p0Expected.*log2(p0Expected));
            normPermEnExpected = permEnExpected/log2(factorial(m));

            out = EN_PermEn(y,m,tau);

            testCase.verifyEqual(out.permEn, permEnExpected, 'AbsTol', 1e-10);
            testCase.verifyEqual(out.normPermEn, normPermEnExpected, 'AbsTol', 1e-10);
        end

        %-------------------------------------------------------------
        % NL_crptool_fnn caching
        %-------------------------------------------------------------
        function test_NL_crptool_fnn_CacheDoesNotCrossContaminate(testCase)
            % NL_crptool_fnn resets the random seed to a fixed value
            % (BF_ResetSeed) before running Marwan's CRPToolbox code, so it is
            % fully deterministic given its inputs, making it safe to cache.
            % Confirm interleaved calls with different series don't get
            % confused with each other's cached values, and that repeated
            % calls give identical (not just similar) results.
            rng(21);
            y1 = randn(500,1);
            y2 = randn(600,1);

            out1a = NL_crptool_fnn(y1);
            out2 = NL_crptool_fnn(y2);
            out1b = NL_crptool_fnn(y1); % should hit the cache

            testCase.verifyEqual(out1a.fnn2, out1b.fnn2, ...
                'A cache hit for y1 must not be contaminated by the intervening call for y2.');
            testCase.verifyEqual(out1a.firstunder05, out1b.firstunder05);
            testCase.verifyNotEqual(out1a.fnn2, out2.fnn2, ...
                'Sanity check: the two different series used here should not coincidentally match exactly.');
        end

        function test_NL_crptool_fnn_CacheEvictionStillCorrect(testCase)
            % Push more distinct series through than the cache holds, then
            % confirm a call for the very first (now-evicted) series still
            % recomputes correctly rather than returning stale/wrong data.
            rng(22);
            N = 500;
            firstSeries = randn(N,1);
            out1a = NL_crptool_fnn(firstSeries);
            for k = 1:6
                NL_crptool_fnn(randn(N,1)); % evict firstSeries from the cache
            end
            out1b = NL_crptool_fnn(firstSeries); % must recompute correctly

            testCase.verifyEqual(out1a.fnn2, out1b.fnn2);
            testCase.verifyEqual(out1a.mdrop, out1b.mdrop);
        end

        %-------------------------------------------------------------
        % DN_OutlierInclude
        %-------------------------------------------------------------
        function test_DN_OutlierInclude_MatchesBruteForceLoop(testCase)
            % DN_OutlierInclude used to loop across the full threshold range
            % and trim off unusable trailing rows (NaN inter-event stats, or
            % less than 2% of data included) afterward; it now breaks out of
            % the loop as soon as either condition is hit (both are
            % monotonic: raising the threshold only shrinks the qualifying
            % point set), and searches only the previous iteration's
            % qualifying indices instead of rescanning the full series each
            % time. Cross-check simple msDt-derived stats against an
            % independent reimplementation of the original brute-force
            % loop + post-hoc trim, for all three thresholdHow modes.
            rng(23);
            N = 1000;
            y = zscore(randn(N,1));
            inc = 0.05; % coarser than the default 0.01, for a faster test

            modes = {'abs','pos','neg'};
            for mi = 1:numel(modes)
                thresholdHow = modes{mi};

                switch thresholdHow
                case 'abs'
                    thr = (0:inc:max(abs(y)));
                    tot = N;
                case 'pos'
                    thr = (0:inc:max(y));
                    tot = sum(y >= 0);
                case 'neg'
                    thr = (0:inc:max(-y));
                    tot = sum(y <= 0);
                end
                msDt = zeros(length(thr),6);
                for i = 1:length(thr)
                    th = thr(i);
                    switch thresholdHow
                    case 'abs'
                        r = find(abs(y) >= th);
                    case 'pos'
                        r = find(y >= th);
                    case 'neg'
                        r = find(y <= -th);
                    end
                    Dt_exc = diff(r);
                    msDt(i,1) = mean(Dt_exc);
                    msDt(i,2) = std(Dt_exc)/sqrt(length(r));
                    msDt(i,3) = length(Dt_exc)/tot*100;
                    msDt(i,4) = median(r)/(N/2)-1;
                    msDt(i,5) = mean(r)/(N/2)-1;
                    msDt(i,6) = std(r)/sqrt(length(r));
                end
                fbi = find(isnan(msDt(:,1)),1,'first');
                if ~isempty(fbi)
                    msDt = msDt(1:fbi-1,:);
                end
                trimthr = 2;
                mj = find(msDt(:,3) > trimthr,1,'last');
                if ~isempty(mj)
                    msDt = msDt(1:mj,:);
                end

                expectedMdtm = mean(msDt(:,1));
                expectedMdrmd = median(msDt(:,4));
                expectedMrstd = std(msDt(:,5));

                out = DN_OutlierInclude(y,thresholdHow,inc);

                testCase.verifyEqual(out.mdtm, expectedMdtm, 'AbsTol', 1e-9, ...
                    sprintf('mdtm mismatch for thresholdHow=''%s''',thresholdHow));
                testCase.verifyEqual(out.mdrmd, expectedMdrmd, 'AbsTol', 1e-9, ...
                    sprintf('mdrmd mismatch for thresholdHow=''%s''',thresholdHow));
                testCase.verifyEqual(out.mrstd, expectedMrstd, 'AbsTol', 1e-9, ...
                    sprintf('mrstd mismatch for thresholdHow=''%s''',thresholdHow));
            end
        end

        %-------------------------------------------------------------
        % BF_CheckToolbox caching
        %-------------------------------------------------------------
        function test_BF_CheckToolbox_CachingDoesNotChangeResult(testCase)
            % BF_CheckToolbox caches the installed-addons list and
            % license('test',...) results (both read-only queries that can't
            % change during a single MATLAB session), since profiling showed
            % installedAddons alone costs ~0.3s per call and this function is
            % called by ~42 operations. Confirm interleaved/repeated calls
            % (which exercise the cache) give results consistent with fresh
            % queries to the underlying Matlab functions, for both the
            % install-check and license-check paths -- without assuming any
            % particular toolbox is absent (environment-independent).
            toolboxes = {'curve_fitting_toolbox','statistics_toolbox','curve_fitting_toolbox'};
            names = {'Curve Fitting Toolbox','Statistics and Machine Learning Toolbox','Curve Fitting Toolbox'};

            for k = 1:numel(toolboxes)
                [outFlag,theName] = BF_CheckToolbox(toolboxes{k},true,true);
                expectedInstalled = any(ismember(matlab.addons.installedAddons().Name, names{k}));
                testCase.verifyEqual(theName, names{k});
                testCase.verifyEqual(outFlag, ~expectedInstalled, ...
                    sprintf('Mismatch for toolbox %s at call %d',toolboxes{k},k));
            end

            % License-check path (doInstallCheck=false):
            [outFlagLic,~] = BF_CheckToolbox('curve_fitting_toolbox',true,false);
            expectedLicensed = license('test','curve_fitting_toolbox');
            testCase.verifyEqual(outFlagLic, ~expectedLicensed);
        end

        %-------------------------------------------------------------
        % CO_TranslateShape
        %-------------------------------------------------------------
        function test_CO_TranslateShape_MatchesExplicitOnesBroadcast(testCase)
            % CO_TranslateShape's 'circle' case used win - ones(2*w+1,1)*ty(...)
            % to broadcast-subtract a row from every row of win; it now relies
            % on implicit broadcasting. Cross-check against the original
            % explicit-ones approach (same fix as CO_Embed2_Shapes).
            rng(25);
            N = 300;
            y = randn(N,1);
            d = 3;

            ty = [(1:N)', y];
            w = floor(d);
            rnge = 1+w:N-w;
            NN = length(rnge);
            npExpected = zeros(NN,1);
            for i = 1:NN
                win = ty(rnge(i)-w:rnge(i)+w,:);
                difwin = win - ones(2*w+1,1)*ty(rnge(i),:);
                npExpected(i) = sum(sum(difwin.^2,2) <= d^2);
            end

            out = CO_TranslateShape(y,'circle',d,'pts');

            testCase.verifyEqual(out.mean, mean(npExpected), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.max, max(npExpected), 'AbsTol', 1e-10);
            testCase.verifyEqual(out.std, std(npExpected), 'AbsTol', 1e-10);
        end

        %-------------------------------------------------------------
        % CO_glscf caching
        %-------------------------------------------------------------
        function test_CO_glscf_CacheDoesNotCrossContaminate(testCase)
            % CO_glscf caches abs(y) (keyed on an exact, NaN-safe match of
            % y), since it's called ~31 times directly (as separate master
            % operations with different alpha/beta/tau on the same series)
            % plus repeatedly from CO_fzcglscf's internal loop. Confirm
            % interleaved calls with different series/parameters don't get
            % confused with each other's cached values, against an
            % independent reference computation.
            rng(26);
            y1 = randn(200,1);
            y2 = randn(250,1);
            alpha = 1; beta = 2; tau = 3;

            referenceGlscf = @(y) refGlscf(y,alpha,beta,tau);

            out1a = CO_glscf(y1,alpha,beta,tau);
            out2 = CO_glscf(y2,alpha,beta,tau);
            out1b = CO_glscf(y1,alpha,beta,tau); % should hit the cache

            testCase.verifyEqual(out1a, referenceGlscf(y1), 'AbsTol', 1e-10);
            testCase.verifyEqual(out2, referenceGlscf(y2), 'AbsTol', 1e-10);
            testCase.verifyEqual(out1b, referenceGlscf(y1), 'AbsTol', 1e-10, ...
                'A cache hit for y1 must not be contaminated by the intervening call for y2.');
        end

        %-------------------------------------------------------------
        % BF_ResetSeed numeric seed support
        %-------------------------------------------------------------
        function test_BF_ResetSeed_NumericSeedIsReproducibleAndDistinct(testCase)
            % BF_ResetSeed originally only accepted 'default'/'none'; a
            % numeric scalar seed was added so surrogate-generating
            % operations can be seeded distinctly from the fixed default.
            BF_ResetSeed(11); a1 = rand(5,1);
            BF_ResetSeed(11); a2 = rand(5,1);
            testCase.verifyEqual(a1, a2, 'A numeric seed should be reproducible.');

            BF_ResetSeed('default'); d1 = rand(5,1);
            testCase.verifyNotEqual(a1, d1, ...
                'A non-zero numeric seed should differ from the ''default'' (seed 0) state.');
        end

        %-------------------------------------------------------------
        % SD_TSTL_surrogates (now native, no TSTOOL dependency)
        %-------------------------------------------------------------
        function test_SD_TSTL_surrogates_DetectsNonlinearity(testCase)
            % SD_TSTL_surrogates used to depend on TSTOOL's signal/tc3/trev;
            % it now generates surrogates via SD_MakeSurrogates and
            % evaluates CO_tc3/CO_trev on each (both already TSTOOL-free
            % reimplementations of the same expressions). Confirm the
            % replacement still does its actual job: a linear (Gaussian
            % noise) series should look unremarkable against its
            % phase-randomized surrogates, while a chaotic (logistic map)
            % series -- whose nonlinear structure phase randomization
            % destroys -- should look highly significant.
            rng(31);
            yLinear = randn(500,1);
            outLinear = SD_TSTL_surrogates(yLinear,1,50,1,'tc3','default');
            testCase.verifyGreaterThan(outLinear.ztestp, 0.05, ...
                'Gaussian noise should not look significantly different from its surrogates.');

            N = 1000;
            yChaotic = zeros(N,1);
            yChaotic(1) = 0.4;
            for i = 2:N
                yChaotic(i) = 3.9*yChaotic(i-1)*(1-yChaotic(i-1));
            end
            outChaotic = SD_TSTL_surrogates(yChaotic,1,50,1,'tc3','default');
            testCase.verifyLessThan(outChaotic.ztestp, 1e-6, ...
                'The logistic map''s nonlinear structure should look highly significant against its surrogates.');
        end

        %-------------------------------------------------------------
        % MF_FitSubsegments 'arma' asymmetric order
        %-------------------------------------------------------------
        function test_MF_FitSubsegments_ArmaAsymmetricOrder(testCase)
            % The 'arma' case's parameter-statistics loops used
            % `for i = 1:order`, where order is a two-element [p,q]
            % vector for this model -- invalid syntax ("Colon operands
            % must be real scalars") whenever order wasn't a scalar,
            % which it never is here. Confirm order(1)/order(2) are now
            % used independently by checking the right number of p_/q_
            % fields are produced for an asymmetric AR/MA order.
            rng(21);
            y = randn(300,1);
            out = MF_FitSubsegments(y,'arma',[2,3],'uniform',[10,0.1]);
            for i = 1:2
                testCase.verifyTrue(isfield(out,sprintf('p_%u_mean',i)));
            end
            testCase.verifyFalse(isfield(out,'p_3_mean'));
            for i = 1:3
                testCase.verifyTrue(isfield(out,sprintf('q_%u_mean',i)));
            end
        end

        %-------------------------------------------------------------
        % NL_TSTL_TakensEstimator (now via TISEAN d2/c2t, no TSTOOL dependency)
        %-------------------------------------------------------------
        function test_NL_TSTL_TakensEstimator_DiscriminatesDimension(testCase)
            % NL_TSTL_TakensEstimator used to depend on TSTOOL's
            % signal/takens_estimator; it now uses TISEAN's d2 and c2t
            % (already relied on by NL_TISEAN_d2.m) instead. Confirm the
            % replacement still does its actual job: high-dimensional
            % (structureless) noise embedded at m should give a dimension
            % estimate close to m, while a logistic map (a classic 1-D
            % chaotic system, known correlation dimension ~1) should give an
            % estimate close to 1, clearly distinguishing it from noise.
            rng(41);
            yNoise = randn(2000,1);
            outNoise = NL_TSTL_TakensEstimator(yNoise,-1,0.05,1,{1,3});
            testCase.verifyGreaterThan(outNoise, 2, ...
                'Structureless noise embedded at m=3 should have a dimension estimate well above 1.');

            N = 2000;
            yChaotic = zeros(N,1);
            yChaotic(1) = 0.4;
            for i = 2:N
                yChaotic(i) = 3.9*yChaotic(i-1)*(1-yChaotic(i-1));
            end
            outChaotic = NL_TSTL_TakensEstimator(yChaotic,-1,0.05,1,{1,3});
            testCase.verifyLessThan(outChaotic, 1.5, ...
                'The logistic map''s known correlation dimension (~1) should give a low estimate.');
            testCase.verifyLessThan(outChaotic, outNoise, ...
                'The logistic map should have a clearly lower dimension estimate than structureless noise.');
        end

        %-------------------------------------------------------------
        % NL_TSTL_GPCorrSum (now via TISEAN d2, no TSTOOL dependency)
        %-------------------------------------------------------------
        function test_NL_TSTL_GPCorrSum_DiscriminatesDimension(testCase)
            % NL_TSTL_GPCorrSum used to depend on TSTOOL's
            % signal/corrsum; it now uses TISEAN's d2 (already relied on
            % by NL_TISEAN_d2.m and NL_TSTL_TakensEstimator.m) instead,
            % fitting its own robust line to d2's raw (r, C(r))
            % correlation-sum output. Confirm the replacement still does
            % its actual job: a logistic map (known correlation dimension
            % ~1) should give a low-slope estimate, clearly below
            % structureless noise's estimate at the same embedding.
            rng(51);
            yNoise = randn(2000,1);
            outNoise = NL_TSTL_GPCorrSum(yNoise,-1,0.5,40,20,{1,3});
            testCase.verifyGreaterThan(outNoise.robfit_a2, 2, ...
                'Structureless noise embedded at m=3 should have a dimension estimate well above 1.');

            N = 2000;
            yChaotic = zeros(N,1);
            yChaotic(1) = 0.4;
            for i = 2:N
                yChaotic(i) = 3.9*yChaotic(i-1)*(1-yChaotic(i-1));
            end
            outChaotic = NL_TSTL_GPCorrSum(yChaotic,-1,0.5,40,20,{1,3});
            testCase.verifyLessThan(outChaotic.robfit_a2, 1.5, ...
                'The logistic map''s known correlation dimension (~1) should give a low estimate.');
            testCase.verifyLessThan(outChaotic.robfit_a2, outNoise.robfit_a2, ...
                'The logistic map should have a clearly lower dimension estimate than structureless noise.');
        end

        function test_NL_TSTL_PoincareSection_DiscriminatesStructure(testCase)
            % NL_TSTL_PoincareSection used to depend on TSTOOL's
            % signal/poincare, which cut a hyperplane orthogonal to the
            % local tangent vector at a chosen reference point. TISEAN's
            % poincare has no such reference-point concept -- it cuts a
            % fixed embedding coordinate at its own mean, in one crossing
            % direction -- so 'ref' ('max'/'min') is now repurposed to pick
            % the crossing direction instead (see the operation's header
            % comment). Confirm the replacement still does its actual job:
            % a smooth chaotic flow (Lorenz) crosses a fixed-threshold
            % hyperplane far less often, per sample, than structureless
            % noise, and its section points cluster into fewer boxes
            % (lower occupancy entropy) than noise's near-uniform fill.
            rng(61);
            yNoise = randn(2000,1);
            outNoise = NL_TSTL_PoincareSection(yNoise,'max',{1,3});

            sigma = 10; rho = 28; beta = 8/3;
            lorenzODE = @(t,v) [sigma*(v(2)-v(1)); v(1)*(rho-v(3))-v(2); v(1)*v(2)-beta*v(3)];
            [~, sol] = ode45(lorenzODE, 0:0.02:400, [1,1,1]);
            yLorenz = sol(2001:end, 1); % discard transient
            yLorenz = zscore(yLorenz(1:5:end)); % downsample
            outLorenz = NL_TSTL_PoincareSection(yLorenz,'max',{1,3});

            testCase.verifyLessThan(outLorenz.pcross, 0.5 * outNoise.pcross, ...
                'The Lorenz flow should cross a fixed-threshold hyperplane far less often than noise.');
            testCase.verifyLessThan(outLorenz.hboxcounts5, outNoise.hboxcounts5, ...
                'The Lorenz section should have lower box-occupancy entropy (more clustered) than noise.');
        end

        function test_NL_TSTL_LargestLyap_DiscriminatesStructure(testCase)
            % NL_TSTL_LargestLyap used to depend on TSTOOL's
            % signal/largelyap; it now uses TISEAN's lyap_r (Rosenstein
            % method) instead, which tracks the same kind of single-
            % nearest-neighbor divergence but has no equivalent of Nref or
            % NNR (always uses all points / exactly one neighbor -- see the
            % operation's header comment). Confirm the replacement still
            % does its actual job: structureless noise has no real
            % correlation between "nearest neighbors", so its divergence
            % curve saturates almost immediately, while a smooth chaotic
            % flow (Lorenz) unfolds its divergence smoothly over many more
            % steps before saturating.
            rng(71);
            yNoise = randn(2000,1);
            outNoise = NL_TSTL_LargestLyap(yNoise,-1,0.1,0.01,3,{1,4});

            sigma = 10; rho = 28; beta = 8/3;
            lorenzODE = @(t,v) [sigma*(v(2)-v(1)); v(1)*(rho-v(3))-v(2); v(1)*v(2)-beta*v(3)];
            [~, sol] = ode45(lorenzODE, 0:0.02:280, [1,1,1]);
            yLorenz = sol(2001:end, 1); % discard transient
            yLorenz = zscore(yLorenz);
            outLorenz = NL_TSTL_LargestLyap(yLorenz,-1,0.1,0.01,3,{1,4});

            testCase.verifyGreaterThan(outLorenz.to09max, 10 * outNoise.to09max, ...
                'The Lorenz flow''s divergence curve should take far longer to saturate than noise''s.');
        end

        function test_NL_TSTL_BoxCorrDim_DiscriminatesStructure(testCase)
            % NL_TSTL_BoxCorrDim used to depend on TSTOOL's signal/corrdim;
            % it now uses TISEAN's boxcount (Renyi entropy of order Q=2.0
            % via box partitioning), taking its per-embedding-dimension
            % entropy increment as the analogue of corrdim's local-dimension
            % matrix (see the operation's header comment). Confirm the
            % replacement still does its actual job using the same
            % noise-vs-logistic-map contrast as the other TSTOOL-derived
            % operations in this file. (A Lorenz-vs-noise contrast was also
            % tried here, but at accessible series lengths the two aren't
            % reliably distinguishable by this particular box-counting
            % statistic -- Lorenz's true correlation dimension, ~2.05, is
            % close enough to embedding dimension m=5 that it hits the same
            % small-epsilon finite-sample sparsity collapse as noise does,
            % at similar length scales; only a much lower-dimensional
            % system like the logistic map (D~1) collapses distinctly
            % later. This is a property of the estimator/sample size, not a
            % migration defect -- the direct old-vs-new TSTOOL comparison
            % during development showed the same qualitative shape.)
            rng(81);
            yNoise = randn(2000,1);
            outNoise = NL_TSTL_BoxCorrDim(yNoise,50,{1,5});

            N = 2000;
            yChaotic = zeros(N,1);
            yChaotic(1) = 0.4;
            for i = 2:N
                yChaotic(i) = 3.9*yChaotic(i-1)*(1-yChaotic(i-1));
            end
            outChaotic = NL_TSTL_BoxCorrDim(yChaotic,50,{1,5});

            testCase.verifyGreaterThan(outChaotic.mediand5, outNoise.mediand5, ...
                ['The logistic map''s deterministic, tightly-clustered embedding should keep boxes ' ...
                 'populated to smaller length scales than noise, giving a higher median entropy increment.']);
        end

        function test_TSTL_localdensity_MatchesBruteForceLoop(testCase)
            % TSTL_localdensity used to depend on TSTOOL's 'localdensity',
            % which the original author noted was "very poorly documented"
            % -- its exact algorithm was never confirmed. It's now a native
            % k-nearest-neighbor local density estimate (density(i) ~
            % 1/r_NNR(i)^m, r_NNR(i) = distance to the NNR-th nearest
            % neighbor outside a Theiler window of "past" samples), using a
            % KD-tree for speed. Since there's no original TSTOOL ground
            % truth to validate against here (unlike the other TSTOOL-
            % derived operations in this file), the right check is that the
            % new native implementation is itself correct: confirm its
            % KD-tree-based fast path matches a brute-force pairwise-
            % distance reference computation with the same Theiler-window
            % exclusion.
            rng(91);
            y = randn(300,1);
            NNR = 4; past = 5;
            out = TSTL_localdensity(y, NNR, past, {1,3});

            Y = BF_Embed(y, 1, 3, 0);
            N_embed = size(Y,1);
            m = size(Y,2);
            locdenBrute = zeros(N_embed,1);
            for i = 1:N_embed
                d = sqrt(sum((Y - Y(i,:)).^2, 2));
                d(abs((1:N_embed)' - i) <= past) = Inf;
                sd = sort(d);
                locdenBrute(i) = 1 / (sd(NNR)^m);
            end

            testCase.verifyEqual(out.meanden, mean(locdenBrute), 'RelTol', 1e-9);
            testCase.verifyEqual(out.stdden, std(locdenBrute), 'RelTol', 1e-9);
            testCase.verifyEqual(out.medianden, median(locdenBrute), 'RelTol', 1e-9);
        end

        %-------------------------------------------------------------
        % Master-operation code string -> function handle equivalence
        %-------------------------------------------------------------
        function test_MasterCodeStringHandleEquivalence(testCase)
            % TS_ComputeMasterLoop used to invoke each master operation's Code
            % string via eval/evalc; it now calls a precompiled function handle
            % (str2func(['@(x,x_z) ' code])) instead, built once in TS_Compute /
            % TS_CalculateFeatureVector. This confirms that switch is behavior-
            % preserving, using real Code strings taken from the mop definition
            % files (FeatureSets/INP_mops_hctsa.txt, INP_mops_catch22.txt).
            rng(6);
            x = randn(300,1);
            x_z = zscore(x);

            sampleCodes = {
                'DN_Mean(x,''norm'')', ...
                'DN_Mean(x,''rms'')', ...
                'DN_Spread(x,''std'')', ...
                'DN_HistogramMode(x_z,5,true,false)', ...
                'catch22_CO_FirstMin_ac(x'')', ...
            };

            for i = 1:numel(sampleCodes)
                code = sampleCodes{i};
                evalResult = eval(code);
                fn = str2func(['@(x,x_z) ',code]);
                handleResult = fn(x,x_z);
                testCase.verifyEqual(handleResult, evalResult, ...
                    sprintf('Handle-based call diverged from eval for code: %s',code));
            end
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

function out = referenceFirstMinAC(y)
    % Independent reference reimplementation of CO_FirstMin's 'ac' branch
    % with default minNotMax=true (the original, uncached algorithm), used to
    % verify CO_FirstMin's cache doesn't cross-contaminate results between
    % different series.
    N = length(y);
    acf_all = CO_AutoCorr(y,[],'Fourier');
    autoCorr = zeros(N-1,1);
    for i = 1:N-1
        autoCorr(i) = acf_all(i+1);
        if isnan(autoCorr(i))
            out = NaN;
            return
        end
        if i==2 && (autoCorr(2) > autoCorr(1))
            out = 1;
            return
        elseif (i > 2) && (autoCorr(i-2) > autoCorr(i-1)) && (autoCorr(i-1) < autoCorr(i))
            out = i-1;
            return
        end
    end
    out = N;
end

function acf = referenceFourierACF(y)
    % Independent reference implementation of CO_AutoCorr's 'Fourier' method
    % (the original, uncached algorithm), used to verify the caching layer
    % added to CO_AutoCorr.m doesn't change its output.
    N = length(y);
    nFFT = 2^(nextpow2(N)+1);
    F = fft(y - mean(y),nFFT);
    F = F.*conj(F);
    acf = ifft(F);
    acf = acf./acf(1);
    acf = real(acf);
    acf = acf(1:N);
end

function glscf = refGlscf(y,alpha,beta,tau)
    % Independent reference implementation of CO_glscf (the original,
    % uncached algorithm: abs(y) recomputed directly, no caching), used to
    % verify CO_glscf's cache doesn't cross-contaminate results between
    % different series.
    y1 = abs(y(1:end-tau));
    y2 = abs(y(1+tau:end));
    glscf = (mean((y1.^alpha).*(y2.^beta)) - mean(y1.^alpha)*mean(y2.^beta)) / ...
                (sqrt(mean(y1.^(2*alpha)) - mean(y1.^alpha)^2) ...
                      * sqrt(mean(y2.^(2*beta)) - mean(y2.^beta)^2));
end
