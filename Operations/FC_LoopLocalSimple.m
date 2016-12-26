function out = FC_LoopLocalSimple(y,forecastMeth)
% FC_LoopLocalSimple    How simple local forecasting depends on window length.
%
% Analyzes the outputs of FC_LocalSimple for a range of local window lengths, l.
% Loops over the length of the data to use for FC_LocalSimple prediction
%
%---INPUTS:
%
% y, the input time series
%
% forecastMeth, the prediction method:
%            (i) 'mean', local mean prediction
%            (ii) 'median', local median prediction
%
%---OUTPUTS:
% Statistics including whether the mean square error increases or decreases,
% testing for peaks, variability, autocorrelation, stationarity, and a fit of
% exponential decay, f(x) = A*exp(Bx) + C, to the variation.

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

% ------------------------------------------------------------------------------
% Check a curve-fitting toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('curve_fitting_toolbox');

% ------------------------------------------------------------------------------
% Also uses xcorr from the signal processing toolbox:
% ------------------------------------------------------------------------------
BF_CheckToolbox('signal_toolbox');

doPlot = 0; % plot outputs to a figure

if nargin < 2 || isempty(forecastMeth)
    forecastMeth = 'mean'; % do mean prediction by default
end

switch forecastMeth
    case 'mean'
        trainLengthRange = (1:10)';

    case 'median'
        trainLengthRange = (1:2:19)';

    otherwise
        error('Unknown prediction method ''%s''',forecastMeth);
end

stats_st = zeros(length(trainLengthRange),5);

for i = 1:length(trainLengthRange)
    switch forecastMeth
        case 'mean'
            outtmp = FC_LocalSimple(y,'mean',trainLengthRange(i));
        case 'median'
            outtmp = FC_LocalSimple(y,'median',trainLengthRange(i));
            % median needs more tweaking
    end
    stats_st(i,1) = outtmp.stderr;
    stats_st(i,2) = outtmp.sws;
    stats_st(i,3) = outtmp.swm;
    stats_st(i,4) = outtmp.ac1;
    stats_st(i,5) = outtmp.ac2;
end

if doPlot
    figure('color','w')
    plot(stats_st,'k')
end

%-------------------------------------------------------------------------------
% Compute statistics from the shapes of the curves
%-------------------------------------------------------------------------------

% (1) root mean square error
% (i) (expect error to decrease with increasing forecast window?:)
out.stderr_chn = mean(diff(stats_st(:,1)))/(range(stats_st(:,1)));
out.stderr_meansgndiff = mean(sign(diff(stats_st(:,1))));

% (ii) is there a peak?
if out.stderr_chn < 1; % on the whole decreasing, as expected
    wigv = max(stats_st(:,1));
    wig = find(stats_st(:,1) == wigv,1,'first');
    if wig ~= 1 && stats_st(wig-1,1) > wigv
        wig = NaN; % maximum is not a local maximum; previous value exceeds it
    elseif wig ~= length(trainLengthRange) && stats_st(wig+1,1) > wigv
        wig = NaN; % maximum is not a local maximum; the next value exceeds it
    end
else
    wigv = min(stats_st(:,1));
    wig = find(stats_st(:,1) == wigv,1,'first');

    if wig ~= 1 && stats_st(wig-1,1) < wigv
        wig = NaN; % maximum is not a local maximum; previous value exceeds it
    elseif wig ~= length(trainLengthRange) && stats_st(wig+1,1) < wigv
        wig = NaN; % maximum is not a local maximum; the next value exceeds it
    end
end
if ~isnan(wig)
    out.stderr_peakpos = wig;
    out.stderr_peaksize = wigv/mean(stats_st(:,1));
else % put NaNs in all the outputs
    out.stderr_peakpos = NaN;
    out.stderr_peaksize = NaN;
end

% (2) Sliding Window Stationarity
out.sws_chn = mean(diff(stats_st(:,2)))/(range(stats_st(:,2)));
out.sws_meansgndiff = mean(sign(diff(stats_st(:,2))));
out.sws_stdn = std(stats_st(:,2))/range(stats_st(:,2));

% Fit exponential decay:
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[range(stats_st(:,2)), -0.5 min(stats_st(:,2))]);
f = fittype('a*exp(b*x)+c','options',s);
[c, gof] = fit(trainLengthRange,stats_st(:,2),f);
out.sws_fexp_a = c.a;
out.sws_fexp_b = c.b; % this is important
out.sws_fexp_c = c.c;
out.sws_fexp_r2 = gof.rsquare; % this is more important!
out.sws_fexp_adjr2 = gof.adjrsquare;
out.sws_fexp_rmse = gof.rmse;

% (3) sliding window mean
out.swm_chn = mean(diff(stats_st(:,3)))/(range(stats_st(:,3)));
out.swm_meansgndiff = mean(sign(diff(stats_st(:,3))));
out.swm_stdn = std(stats_st(:,3))/range(stats_st(:,3));

% (4) AC1
out.ac1_chn = mean(diff(stats_st(:,4)))/(range(stats_st(:,4)));
out.ac1_meansgndiff = mean(sign(diff(stats_st(:,4))));
out.ac1_stdn = std(stats_st(:,4))/range(stats_st(:,4));

% (5) AC2
out.ac2_chn = mean(diff(stats_st(:,5)))/(range(stats_st(:,5)));
out.ac2_meansgndiff = mean(sign(diff(stats_st(:,5))));
out.ac2_stdn = std(stats_st(:,5))/range(stats_st(:,5));

end
