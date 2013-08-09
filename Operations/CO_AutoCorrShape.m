% CO_AutoCorrShape
% 
% Outputs a set of statistics summarizing how the autocorrelation function
% changes with the time lag, tau.
% Outputs include the number of peaks, and autocorrelation in the
% autocorrelation function itself.
% 
% INPUTS:
% y, the input time series
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = CO_AutoCorrShape(y)
% Ben Fulcher, 2009

doplot = 0; % plot outputs from this function
N = length(y); % length of the time series

% Only look up to when two consecutive values are under the threshold for
% significance:
th = 2/sqrt(N); % significance threshold, th

% Initialize autocorrelation function
acf = zeros(N,1);

% Calculate the autocorrelation function, up to a maximum lag, length of
% time series (hopefully it's cropped by then)
for i = 1:N
    acf(i) = CO_AutoCorr(y,i-1); % *** NOTE THIS! *** acf vector indicies are not lags
    if i > 2 && abs(acf(i)) < th && abs(acf(i-1)) < th
       acf = acf(1:i-2);
       break
    end
end

Nac = length(acf);

out.Nac = length(acf); % the distance the acf lasts until significance is 'drowned out' (by my definition)

% Count peaks
dacf = diff(acf);
ddacf = diff(dacf);
extrr = BF_sgnchange(dacf,1);
sdsp = ddacf(extrr);
maxr = extrr(sdsp < 0);
minr = extrr(sdsp > 0);
nmaxr = length(maxr);
nminr = length(minr);

% Number of local minima
out.nminima = sum(sdsp > 0);
out.meanminima = mean(sdsp(sdsp > 0));

% Proportion of local maxima
out.nmaxima = sum(sdsp < 0);
out.meanmaxima = abs(mean(sdsp(sdsp < 0))); % must be negative: make it positive

% Proportion of extrema
out.nextrema = length(sdsp);
out.pextrema = length(sdsp)/Nac;

% Mean of the ACF
out.meanacf = mean(acf);
out.meanabsacf = mean(abs(acf));

% Correlations between extrema
if nmaxr > 4 % need at least 5 points to do this
    out.maximaspread = std(diff(maxr)); % spread of inter-maxima intervals
    out.ac1maxima = CO_AutoCorr(acf(maxr),1);
else % less than 5 points, return NaNs:
    out.maximaspread = NaN;
    out.ac1maxima = NaN;
end
if nminr > 4 % need at least 5 points to do this
    out.minimaspread = std(diff(minr)); % spread of inter-minima intervals
    out.ac1minima = CO_AutoCorr(acf(minr),1);
else % less than 5 points, return NaNs:
    out.minimaspread = NaN;
    out.ac1minima = NaN;
end

% Autocorrelation of the ACF
out.ac1 = CO_AutoCorr(acf,1);
out.ac2 = CO_AutoCorr(acf,2);
out.ac3 = CO_AutoCorr(acf,3);
out.actau = CO_AutoCorr(acf,CO_FirstZero(acf,'ac'));


if Nac > 3 % Need at least four points to fit exponential
    
    %% Fit exponential decay to absolute ACF:
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1, -0.5]);
    f = fittype('a*exp(b*x)','options',s);
    b = 1;
    try
        [c, gof] = fit((1:Nac)',abs(acf),f);
    catch
        b = 0;
    end
    if b == 1
        out.fexpabsacf_a = c.a;
        out.fexpabsacf_b = c.b; % this is important
        out.fexpabsacf_r2 = gof.rsquare; % this is more important!
        out.fexpabsacf_adjr2 = gof.adjrsquare;
        out.fexpabsacf_rmse = gof.rmse;
    
        expfit = c.a*exp(c.b*[1:Nac]');
        res = abs(acf)-expfit;
        out.fexpabsacf_varres = var(res);
    else % fit failed -- return NaNs
        out.fexpabsacf_a = NaN;
        out.fexpabsacf_b = NaN;
        out.fexpabsacf_r2 = NaN;
        out.fexpabsacf_adjr2 = NaN;
        out.fexpabsacf_rmse = NaN;
        out.fexpabsacf_varres = NaN;
    end
    
    %% fit linear to local maxima
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[-0.1 1]);
    f = fittype('a*x+b','options',s);
    if doplot
        figure('color','w');
        plot(maxr,acf(maxr),'ok');
    end
    
    b = 1;
    try [c, gof] = fit(maxr,acf(maxr),f);
    catch
        b = 0;
    end
    if b == 1; % Fit was successful
        out.flinlmxacf_a = c.a;
        out.flinlmxacf_b = c.b;
        out.flinlmxacf_r2 = gof.rsquare;
        out.flinlmxacf_adjr2 = gof.adjrsquare;
        out.flinlmxacf_rmse = gof.rmse;
    else % Fit failed -- return NaNs
        out.flinlmxacf_a = NaN;
        out.flinlmxacf_b = NaN;
        out.flinlmxacf_r2 = NaN;
        out.flinlmxacf_adjr2 = NaN;
        out.flinlmxacf_rmse = NaN;
    end
else
    out.fexpabsacf_a = NaN;
    out.fexpabsacf_b = NaN;
    out.fexpabsacf_r2 = NaN;
    out.fexpabsacf_adjr2 = NaN;
    out.fexpabsacf_rmse = NaN;
    out.fexpabsacf_varres = NaN;
    out.flinlmxacf_a = NaN;
    out.flinlmxacf_b = NaN;
    out.flinlmxacf_r2 = NaN;
    out.flinlmxacf_adjr2 = NaN;
    out.flinlmxacf_rmse = NaN;
end


end