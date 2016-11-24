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
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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
    % Stop at double the first time ACF hits within a threshold of ~zero
    stopWhen = 10;
end

% ------------------------------------------------------------------------------
% Check a curve-fitting toolbox license is available
% ------------------------------------------------------------------------------
BF_CheckToolbox('curve_fitting_toolbox');

doPlot = 0; % plot outputs from this function
N = length(y); % length of the time series

% Only look up to when two consecutive values are under the significance threshold
th = 2/sqrt(N); % significance threshold, th

%-------------------------------------------------------------------------------
% Calculate the autocorrelation function, up to a maximum lag, length of
% time series (hopefully it's cropped by then)
%-------------------------------------------------------------------------------
% Compute the autocorrelation function, acf
acf = zeros(N,1);

if isnumeric(stopWhen)
    acf = CO_AutoCorr(y,0:stopWhen,'Fourier');
    Ndrown = stopWhen;
elseif ischar(stopWhen) % compute ACF up to a given threshold:
    switch stopWhen
    case 'drown'
        % Stop when ACF ~ 0
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
        Ndrown = 0; % the point at which ACF ~ 0
        for i = 1:N
            acf(i) = CO_AutoCorr(y,i-1,'Fourier'); % *** NOTE acf vector indicies are not lags
            if (Ndrown > 0) && (i==Ndrown*2)
                acf = acf(1:i);
                break
            elseif (i > 1) && (abs(acf(i)) < th)
                Ndrown = i;
            end
        end
    end
end

if doPlot
    f = figure('color','w'); hold on
    plot(acf,'o-k')
    plot([1,length(acf)],th*ones(2,1),'--k')
    plot([1,length(acf)],-th*ones(2,1),'--k')
    xlabel('time delay')
end

out.Nac = Ndrown; % the distance the acf lasts until significance is 'drowned out' (by my definition)

Nac = length(acf);

%-------------------------------------------------------------------------------
% Basic stats on the ACF
%-------------------------------------------------------------------------------
out.meanacf = mean(acf);
out.meanabsacf = mean(abs(acf));

% Autocorrelation of the ACF
out.ac1 = CO_AutoCorr(acf,1,'Fourier');
out.actau = CO_AutoCorr(acf,CO_FirstZero(acf,'ac'),'Fourier');

%-------------------------------------------------------------------------------
% Local extrema
%-------------------------------------------------------------------------------
dacf = diff(acf);
ddacf = diff(dacf);
extrr = BF_sgnchange(dacf,1);
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

if Nac > 3 % Need at least four points to fit exponential

    %-------------------------------------------------------------------------------
    %% Fit exponential decay to absolute ACF:
    %-------------------------------------------------------------------------------
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1, -0.5]);
    f = fittype('a*exp(b*x)','options',s);
    fitSuccess = 0;
    try
        [c, gof] = fit((1:Nac)',abs(acf),f); fitSuccess = 1;
    end
    if fitSuccess % Fit was successful
        out.fexpabsacf_a = c.a;
        out.fexpabsacf_b = c.b; % this is important
        out.fexpabsacf_r2 = gof.rsquare; % this is more important!
        out.fexpabsacf_adjr2 = gof.adjrsquare;
        out.fexpabsacf_rmse = gof.rmse;

        expfit = c.a*exp(c.b*(1:Nac)');
        res = abs(acf)-expfit;
        out.fexpabsacf_stdres = std(res);

    else % fit failed -- return NaNs
        out.fexpabsacf_a = NaN;
        out.fexpabsacf_b = NaN;
        out.fexpabsacf_r2 = NaN;
        out.fexpabsacf_adjr2 = NaN;
        out.fexpabsacf_rmse = NaN;
        out.fexpabsacf_stdres = NaN;
    end

    % %-------------------------------------------------------------------------------
    % %% Fit linear to local maxima
    % %-------------------------------------------------------------------------------
    % s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[-0.1 1]);
    % f = fittype('a*x+b','options',s);
    % if doPlot
    %     figure('color','w');
    %     plot(maxr,acf(maxr),'ok');
    % end
    %
    % fitSuccess = 0;
    % try [c, gof] = fit(maxr,acf(maxr),f); fitSuccess = 1;
    % end
    % if fitSuccess % Fit was successful
    %     out.flinlmxacf_a = c.a;
    %     out.flinlmxacf_b = c.b;
    %     out.flinlmxacf_r2 = gof.rsquare;
    %     out.flinlmxacf_adjr2 = gof.adjrsquare;
    %     out.flinlmxacf_rmse = gof.rmse;
    % else % Fit failed -- return NaNs
    %     out.flinlmxacf_a = NaN;
    %     out.flinlmxacf_b = NaN;
    %     out.flinlmxacf_r2 = NaN;
    %     out.flinlmxacf_adjr2 = NaN;
    %     out.flinlmxacf_rmse = NaN;
    % end
else
    out.fexpabsacf_a = NaN;
    out.fexpabsacf_b = NaN;
    out.fexpabsacf_r2 = NaN;
    out.fexpabsacf_adjr2 = NaN;
    out.fexpabsacf_rmse = NaN;
    out.fexpabsacf_stdres = NaN;
    % out.flinlmxacf_a = NaN;
    % out.flinlmxacf_b = NaN;
    % out.flinlmxacf_r2 = NaN;
    % out.flinlmxacf_adjr2 = NaN;
    % out.flinlmxacf_rmse = NaN;
end


end
