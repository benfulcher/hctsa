function varargout = BF_GARCHCache(action, varargin)
% BF_GARCHCache   Work shared by MF_GARCHfit and MF_GARCHcompare on the same series.
%
% The operations fit GARCH models to the same preprocessed series (the same input, the
% same preprocessing) and repeat the preprocessing, the data construction, the tests of
% the series and the Gaussian GARCH(P,Q) fit. Cache the last few results, keyed on the
% exact inputs (all deterministic: the random numbers of the preprocessing come from the
% portable generator with a fixed seed; randomSeed 'none' is never cached).
%
% [y, D] = BF_GARCHCache('prep', x, preProc, randomSeed)
%       the whitened, z-scored series y and its data structure D (see BF_GARCHFit)
% [pEngle, pLbq] = BF_GARCHCache('tests', y)
%       Engle's ARCH test and the Ljung-Box Q-test of the squares, lags 1:20
% fit = BF_GARCHCache('gauss', y, D, P, Q, isGJR)
%       the Gaussian maximum-likelihood fit of BF_GARCHFit (grid starts only)

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

persistent prepX prepKey prepY prepD fitY fitTests fitGauss
maxCacheEntries = 3;
if isempty(prepX)
    prepX = {}; prepKey = {}; prepY = {}; prepD = {};
    fitY = {}; fitTests = {}; fitGauss = {};
end

switch action
case 'prep'
    [x, preProc, randomSeed] = deal(varargin{:});
    cacheable = ~(ischar(randomSeed) && strcmp(randomSeed, 'none'));
    if cacheable
        if strcmp(preProc, 'ar')
            key = {preProc, BF_RandomSeed(randomSeed)}; % (the seed only matters here)
        else
            key = {preProc, []};
        end
        for ci = 1:numel(prepX)
            if isequal(key, prepKey{ci}) && isequaln(x, prepX{ci})
                varargout = {prepY{ci}, prepD{ci}};
                return
            end
        end
    end
    y = BF_Whiten(x, preProc, false, randomSeed);
    y = zscore(y); % z-score the time series (after whitening)
    D = BF_GARCHFit('data', y);
    if cacheable
        prepX = [{x}, prepX]; prepKey = [{key}, prepKey];
        prepY = [{y}, prepY]; prepD = [{D}, prepD];
        n = min(maxCacheEntries, numel(prepX));
        prepX = prepX(1:n); prepKey = prepKey(1:n); prepY = prepY(1:n); prepD = prepD(1:n);
    end
    varargout = {y, D};
case 'tests'
    y = varargin{1};
    ci = findY(y);
    if ci == 0 || isempty(fitTests{ci})
        pEngle = BF_ARCHTest(y, 20);
        pLbq = BF_LjungBoxTest(y.^2, 20);
        if ci == 0, ci = addY(y); end
        fitTests{ci} = {pEngle, pLbq};
    end
    varargout = fitTests{ci};
case 'gauss'
    [y, D, P, Q, isGJR] = deal(varargin{:});
    ci = findY(y);
    key = [P, Q, double(logical(isGJR))];
    if ci > 0
        for k = 1:size(fitGauss{ci}, 1)
            if isequal(fitGauss{ci}{k, 1}, key)
                varargout{1} = fitGauss{ci}{k, 2};
                return
            end
        end
    end
    fit = BF_GARCHFit('mle', D, P, Q, isGJR, false);
    if ci == 0, ci = addY(y); end
    fitGauss{ci} = [fitGauss{ci}; {key, fit}];
    varargout{1} = fit;
otherwise
    error('Unknown action ''%s''', action);
end

    function ci = findY(y)
        ci = 0;
        for c = 1:numel(fitY)
            if isequaln(y, fitY{c}), ci = c; return, end
        end
    end

    function ci = addY(y)
        fitY = [{y}, fitY]; fitTests = [{[]}, fitTests]; fitGauss = [{cell(0, 2)}, fitGauss];
        if numel(fitY) > maxCacheEntries
            fitY(end) = []; fitTests(end) = []; fitGauss(end) = [];
        end
        ci = 1;
    end

end
