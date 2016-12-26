function out = PH_Walker(y,walkerRule,walkerParams)
% PH_Walker Simulates a hypothetical walker moving through the time domain.
%
% The hypothetical particle (or 'walker') moves in response to values of the
% time series at each point.
%
% Outputs from this operation are summaries of the walker's motion, and
% comparisons of it to the original time series.
%
%---INPUTS:
%
% y, the input time series
%
% walkerRule, the kinematic rule by which the walker moves in response to the
%             time series over time:
%
%            (i) 'prop': the walker narrows the gap between its value and that
%                        of the time series by a given proportion p.
%                        walkerParams = p;
%
%            (ii) 'biasprop': the walker is biased to move more in one
%                         direction; when it is being pushed up by the time
%                         series, it narrows the gap by a proportion p_{up},
%                         and when it is being pushed down by the time series,
%                         it narrows the gap by a (potentially different)
%                         proportion p_{down}. walkerParams = [pup,pdown].
%
%            (iii) 'momentum': the walker moves as if it has mass m and inertia
%                         from the previous time step and the time series acts
%                         as a force altering its motion in a classical
%                         Newtonian dynamics framework. [walkerParams = m], the mass.
%
%             (iv) 'runningvar': the walker moves with inertia as above, but
%                         its values are also adjusted so as to match the local
%                         variance of time series by a multiplicative factor.
%                         walkerParams = [m,wl], where m is the inertial mass and wl
%                         is the window length.
%
% walkerParams, the parameters for the specified walkerRule, explained above.
%
%---OUTPUTS: include the mean, spread, maximum, minimum, and autocorrelation of
% the walker's trajectory, the number of crossings between the walker and the
% original time series, the ratio or difference of some basic summary statistics
% between the original time series and the walker, an Ansari-Bradley test
% comparing the distributions of the walker and original time series, and
% various statistics summarizing properties of the residuals between the
% walker's trajectory and the original time series.

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
%% Preliminaries
% ------------------------------------------------------------------------------
doPlot = 0; % plot outputs to figure
N = length(y); % the length of the input time series, y

% ------------------------------------------------------------------------------
%% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(walkerRule)
    walkerRule = 'prop'; % default
end

if nargin < 3 || isempty(walkerParams)
    % Set default parameter for this type of walker dynamics
    switch walkerRule
    case 'prop'
        walkerParams = 0.5;
    case 'biasprop'
        walkerParams = [0.1, 0.2];
    case 'momentum'
        walkerParams = 2;
    case 'runningvar'
        walkerParams = [1.5, 50];
    end
end

% ------------------------------------------------------------------------------
%% (1) Walk
% ------------------------------------------------------------------------------

w = zeros(N,1); % the walker's trajectory, w

switch walkerRule
    case 'prop'
        % walker starts at zero and narrows the gap between its position
        % and the time series value at that point by the proportion given
        % in walkerParams, to give the value at the subsequent time step
        p = walkerParams;

        w(1) = 0; % start at zero
        for i = 2:N
            w(i) = w(i-1) + p*(y(i-1)-w(i-1));
        end

    case 'biasprop'
        % walker is biased in one or the other direction (i.e., prefers to
        % go up, or down). Requires a vector of inputs: [p_up, p_down]
        pup = walkerParams(1);
        pdown = walkerParams(2);

        w(1) = 0;
        for i = 2:N
            if y(i) > y(i-1) % time series increases
                w(i) = w(i-1) + pup*(y(i-1)-w(i-1));
            else
                w(i) = w(i-1) + pdown*(y(i-1)-w(i-1));
            end
        end

    case 'momentum'
        % walker moves as if it had inertia from the previous time step,
        % i.e., it 'wants' to move the same amount; the time series acts as
        % a force changing its motion
        m = walkerParams(1); % 'inertial mass'
%         F=walkerParams(2); % weight of 'force' from time series

        w(1) = y(1);
        w(2) = y(2);
        for i = 3:N
            w_inert = w(i-1) + (w(i-1)-w(i-2));
%             w(i)=w_inert+(y(i-1)-w(i-1))/m; % dissipative term
            w(i) = w_inert + (y(i)-w_inert)/m; % dissipative term
            % equation of motion (s-s_0=ut+F/m*t^2)
            % where the 'force' F is the change in the original time series
            % at that point
        end

    case 'runningvar'
        % walker moves with momentum defined by amplitude of past values in
        % a given length window
        m = walkerParams(1); % 'inertial mass'
        wl = walkerParams(2); % window length

        w(1) = y(1);
        w(2) = y(2);
        for i = 3:N
            w_inert = w(i-1) + (w(i-1)-w(i-2));
            w_mom = w_inert + (y(i)-w_inert)/m; % dissipative term from time series
            if i > wl
                w(i) = w_mom*(std(y(i-wl:i))/std(w(i-wl:i))); % adjust by local standard deviation
            else
                w(i) = w_mom;
            end
        end

    otherwise
        error('Unknown method ''%s'' for simulating walker on the time series', walkerRule)
end

% ------------------------------------------------------------------------------
%% PLOT WALKER AND ORIGINAL TIME SERIES TOGETHER:
% ------------------------------------------------------------------------------
if doPlot
    lw = 1; % set the line width for plotting
    figure('color','w'); box('on'); hold on;
    c = BF_getcmap('set1',3,1);
    plot(y,'.-k','LineWidth',lw); % original time series
    plot(w,'.-','color',c{1},'LineWidth',lw); % walker
    plot([1,length(w)],ones(2,1)*mean(w),'color',c{2},'LineWidth',2); % mean
    % running variance:
    stds = nan(N,2);
    for i = wl+1:N
        stds(i,1) = std(y(i-wl:i));
        stds(i,2) = std(w(i-wl:i));
    end
    % plot(stds(:,1),':r'); % this is the time series
    plot(stds(:,1)./stds(:,2),'color',c{3},'LineWidth',lw); % this is the adjustment factor
    % means = zeros(N,1);
    % for i = 1:N
    %     means(i) = mean(w(1:i));
    % end
    % plot(means,'g')
    % plot(y-w,'m'); % residual
    legend('y','walker','mean: walker','localvariancefactor','accumulative walker mean')
end

% ------------------------------------------------------------------------------
%% (2) Statistics on the walk
% ------------------------------------------------------------------------------
% (i) The walk itself
out.w_mean = mean(w);
out.w_median = median(w);
out.w_std = std(w);
out.w_ac1 = CO_AutoCorr(w,1,'Fourier');
out.w_ac2 = CO_AutoCorr(w,2,'Fourier');
out.w_tau = CO_FirstZero(w,'ac');
out.w_min = min(w);
out.w_max = max(w);
out.w_propzcross = sum(w(1:end-1).*w(2:end) < 0) / (N-1);
% fraction of time series length that walker crosses time series

% (ii) Differences between the walk at signal
out.sw_meanabsdiff = mean(abs(y-w));
out.sw_taudiff = CO_FirstZero(y,'ac') - CO_FirstZero(w,'ac');
out.sw_stdrat = std(w)/std(y); % will be the same as w_std for z-scored signal
out.sw_ac1rat = out.w_ac1/CO_AutoCorr(y,1);
out.sw_minrat = min(w)/min(y);
out.sw_maxrat = max(w)/max(y);
out.sw_propcross = sum((w(1:end-1)-y(1:end-1)).*(w(2:end)-y(2:end)) < 0)/(N-1);
% fraction of time series length that walker crosses time series

% test from same distribution: Ansari-Bradley test
% (this may not be valid given the dependence of w on y, and the
% properties of the null hypothesis itself... But this is the name of the game!)
[h, pval, stats] = ansaribradley(w,y);
out.sw_ansarib_pval = pval; % p-value from the test
% out.sw_ansarib_W = stats.W; % W (test statistic)
% out.sw_ansarib_Wstar = stats.Wstar; % Approximate normal statistic
% test statistics are length dependent. Remove.

r = linspace(min(min(y),min(w)),max(max(y),max(w)),200); % make range of ksdensity uniform across all subsegments
dy = ksdensity(y,r); dw = ksdensity(w,r); % the kernel-smoothed distributions
out.sw_distdiff = sum(abs(dy-dw));

% (iii) Looking at residuals between time series and walker
res = w - y;
[h, pval] = runstest(res); % runs test
out.res_runstest = pval;
out.res_swss5_1 = SY_SlidingWindow(res,'std','std',5,1); % sliding window stationarity
out.res_ac1 = CO_AutoCorr(res,1); % auto correlation at lag-1

end
