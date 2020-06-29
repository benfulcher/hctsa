function out = CO_AutoCorrShape(y,stopWhen)
% CO_AutoCorrShape   How the autocorrelation function changes with the time lag.
%
% Outputs include the number of peaks, and autocorrelation in the
% autocorrelation function (ACF) itself.
%
%---INPUTS:
% y, the input time series
% stopWhen, the criterion for the maximum lag to measure the ACF up to.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
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

%-------------------------------------------------------------------------------
% INPUTS:
%-------------------------------------------------------------------------------
if nargin < 2
    % Stop looking at a given time lag
    stopWhen = 'posDrown';
end

% ------------------------------------------------------------------------------
% Check a curve-fitting toolbox license is available
% ------------------------------------------------------------------------------
BF_CheckToolbox('curve_fitting_toolbox');

doPlot = false; % plot outputs from this function
N = length(y); % length of the time series

% Only look up to when two consecutive values are under the significance threshold
th = 2/sqrt(N); % significance threshold, th

%-------------------------------------------------------------------------------
% Calculate the autocorrelation function, up to a maximum lag, length of
% time series (hopefully it's cropped by then)
%-------------------------------------------------------------------------------
% Compute the autocorrelation function, acf
acf = zeros(N,1);

% At what lag does the acf drop to zero, Ndrown (by my definition)?
if isnumeric(stopWhen)
    acf = CO_AutoCorr(y,0:stopWhen,'Fourier');
    Ndrown = stopWhen;
elseif ischar(stopWhen) % compute ACF up to a given threshold:
    Ndrown = 0; % the point at which ACF ~ 0
    switch stopWhen
    case 'posDrown'
        % Stop when ACF drops below threshold, th
        for i = 1:N
            acf(i) = CO_AutoCorr(y,i-1,'Fourier'); % *** NOTE THIS! *** acf vector indicies are not lags
            if isnan(acf(i))
                warning('Weird time series (constant?)');
                out = NaN; return
            end
            if acf(i) < th
                % Ensure ACF is all positive:
                if acf(i) > 0
                    Ndrown = i;
                    acf = acf(1:i);
                else
                    Ndrown = i-1;
                    acf = acf(1:i-1);
                end
                break
            end
        end
        % This should yield the initial, positive portion of the ACF
        assert(all(acf > 0));
    case 'drown'
        % Stop when ACF is very close to 0 (within threshold, th = 2/sqrt(N))
        for i = 1:N
            acf(i) = CO_AutoCorr(y,i-1,'Fourier'); % *** NOTE THIS! *** acf vector indicies are not lags
            if (i > 1) && (abs(acf(i)) < th)
                Ndrown = i;
                acf = acf(1:i);
                break
            end
        end
    case 'doubleDrown'
        % Stop at 2*tau, where tau is the lag where ACF ~ 0 (within 1/sqrt(N) threshold)
        for i = 1:N
            acf(i) = CO_AutoCorr(y,i-1,'Fourier'); % *** NOTE acf vector indicies are not lags
            if (Ndrown > 0) && (i==Ndrown*2)
                acf = acf(1:i);
                break
            elseif (i > 1) && (abs(acf(i)) < th)
                Ndrown = i;
            end
        end
    otherwise
        error('Unknown ACF decay criterion: ''%s''',stopWhen);
    end
end

% Check for good behavior:
if any(isnan(acf))
    % This is an anomalous time series (e.g., all constant, or conatining NaNs)
    out = NaN;
    return
end

if doPlot
    f = figure('color','w'); hold on
    plot(acf,'o-k')
    plot([1,length(acf)],th*ones(2,1),'--k')
    plot([1,length(acf)],-th*ones(2,1),'--k')
    xlabel('time delay')
end

out.Nac = Ndrown;

Nac = length(acf);

%-------------------------------------------------------------------------------
% Basic stats on the ACF
%-------------------------------------------------------------------------------
out.sumacf = sum(acf);
out.meanacf = mean(acf);
if ~strcmp(stopWhen,'posDrown')
    % Can have negative entries:
    out.meanabsacf = mean(abs(acf));
    out.sumabsacf = sum(abs(acf));
end

% Autocorrelation of the ACF
minPointsForACFofACF = 5; % can't take lots of complex stats with fewer than this
if Nac > minPointsForACFofACF
    out.ac1 = CO_AutoCorr(acf,1,'Fourier');
    if all(acf > 0)
        out.actau = NaN;
    else
        out.actau = CO_AutoCorr(acf,CO_FirstZero(acf,'ac'),'Fourier');
    end
else
    out.ac1 = NaN;
    out.actau = NaN;
end

%-------------------------------------------------------------------------------
% Local extrema
%-------------------------------------------------------------------------------
dacf = diff(acf);
ddacf = diff(dacf);
extrr = BF_SignChange(dacf,1);
sdsp = ddacf(extrr);
% maxr = extrr(sdsp < 0);
% minr = extrr(sdsp > 0);
% nmaxr = length(maxr);
% nminr = length(minr);

% Proportion of local minima
out.nminima = sum(sdsp > 0);
out.meanminima = mean(sdsp(sdsp > 0));

% Proportion of local maxima
out.nmaxima = sum(sdsp < 0);
out.meanmaxima = abs(mean(sdsp(sdsp < 0))); % must be negative: make it positive

% Proportion of extrema
out.nextrema = length(sdsp);
out.pextrema = length(sdsp)/Nac;

% Correlations between extrema
% if nmaxr > 4 % need at least 5 points to do this
%     out.maximaspread = std(diff(maxr)); % spread of inter-maxima intervals
%     out.ac1maxima = CO_AutoCorr(acf(maxr),1,'Fourier');
% else % less than 5 points, return NaNs:
%     out.maximaspread = NaN;
%     out.ac1maxima = NaN;
% end
% if nminr > 4 % need at least 5 points to do this
%     out.minimaspread = std(diff(minr)); % spread of inter-minima intervals
%     out.ac1minima = CO_AutoCorr(acf(minr),1,'Fourier');
% else % less than 5 points, return NaNs:
%     out.minimaspread = NaN;
%     out.ac1minima = NaN;
% end

%-------------------------------------------------------------------------------
% FIT EXPONENTIAL DECAY (only for 'posDrown', and if there are enough points)
%-------------------------------------------------------------------------------
% Should probably only do this up to the first zero crossing...
fitSuccess = false;
minPointsToFitExp = 4; % (need at least four points to fit exponential)
if strcmp(stopWhen,'posDrown') & (Nac >= minPointsToFitExp)
    %-------------------------------------------------------------------------------
    %% Fit exponential decay to (absolute) ACF:
    % (kind of only makes sense for the first positive period)
    %-------------------------------------------------------------------------------
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',0.5);
    f = fittype('exp(-b*x)','options',s);
    try
        [c, gof] = fit((0:Nac-1)',acf,f);
        fitSuccess = true;
    end
end
if fitSuccess % Fit was successful
    out.decayTimescale = 1./c.b; % this is important
    out.fexpacf_r2 = gof.rsquare; % this is more important!
    % out.fexpacf_adjr2 = gof.adjrsquare;
    % out.fexpacf_rmse = gof.rmse;

    expfit = exp(c.b*(0:Nac-1)');
    residuals = acf - expfit;
    out.fexpacf_stdres = std(residuals);

else % fit inappropriate (or failed): return NaNs for the relevant stats
    out.decayTimescale = NaN;
    out.fexpacf_r2 = NaN;
    % out.fexpacf_adjr2 = NaN;
    % out.fexpacf_rmse = NaN;
    out.fexpacf_stdres = NaN;
end


end
