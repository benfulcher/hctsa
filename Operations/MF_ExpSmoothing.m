function out = MF_ExpSmoothing(x,ntrain,alpha)
% MF_ExpSmoothing   Exponential smoothing time-series prediction model.
%
% Fits an exponential smoothing model to the time series using a training set to
% fit the optimal smoothing parameter, alpha, and then applies the result to the
% try to predict the rest of the time series.
%
% cf. "The Analysis of Time Series", C. Chatfield, CRC Press LLC (2004)
%
% Code is adapted from original code provided by Siddharth Arora:
% Siddharth.Arora@sbs.ox.ac.uk
%
%---INPUTS:
% x, the input time series
%
% ntrain, the number of samples to use for training (can be a proportion of the
%           time-series length)
%
% alpha, the exponential smoothing parameter
%
%---OUTPUTS: include the fitted alpha, and statistics on the residuals from the
% prediction phase.
%
% Future alteration could take a number of training sets and average to some
% optimal alpha, for example, rather than just fitting it in an initial portion
% of the time series.

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

doPlot = 0; % plot outputs
N = length(x); % the length of the time series

%-------------------------------------------------------------------------------
% Check Inputs:
%-------------------------------------------------------------------------------

% (*) ntrain -- either the number or proportion of training points
if nargin < 2 || isempty(ntrain)
    ntrain = min(100,N); % if the time series is shorter than 100 samples(!)
end
if (ntrain > 0) && (ntrain < 1) % given training set length as proportion
                            % of the time series length
    ntrain = floor(N*ntrain);
end

% Check training set sizes:
mintrain = 100; % Minimum training set size
maxtrain = 1000; % Maximum training set size

if ntrain > maxtrain; % larger than maximum training set size
    fprintf(1,'Training set size reduced from %u to maximum of 1000 samples.\n',ntrain);
    ntrain = 1000;
end
if ntrain < mintrain; % smaller than minimum training set size
    fprintf(1,'Training set size increased from %u to minimum of 100.\n',ntrain);
    ntrain = 100;
end

if N < ntrain % time series shorter than the size of the training set
    fprintf(1,'Time Series too short for exponential smoothing\n')
    out = NaN; return
end

% (*) alpha, the smoothing parameter
if nargin < 3 || isempty(alpha)
    alpha = 'best';
end


% Exponential smoothing
% Using: S(t+1) = A.X(t+1) + (1-A).S(t) where S(1) = X(1)
% Finding optimum parameter A [0,1] using RMSE

% n1 = length(x);
% count = 1;

if strcmp(alpha,'best')
    %% (*) Optimize alpha (*)
    % optimize alpha over the training set xtrain. This is a choice to only use
    % the first section of the time series. The length of the training set is
    % set in the input
    % optimize alpha based on the defined training set of time series values

    % Training set xtrain
    xtrain = x(1:ntrain);

    bruteforce = 0;

    if bruteforce
        % set a range of alpha and find minimum
        alphar = (0.1:0.1:1); % the exponential smoothing parameter
        nalpha = length(alphar);
        rmses = zeros(nalpha,1);

        for k = 1:nalpha;
            a = alphar(k);

            % Loop for rolling window
            xf = SUB_fit_exp_smooth(xtrain,a);

            % Issue forecasts
            fore = xf(3:end);
            orig = xtrain(3:end);
            rmses(j) = sqrt(mean((fore - orig).^2)); % compute rmse
%             rmses(k) = rmse(fore,orig);
%             rmse_n(k,2) = a;

            % plot(orig);hold on;plot(fore,'r');hold off
            input(num2str(a));

        %     count = count + 1;
            clear xf;
            clear fore;

            % ++BF -- halts unnecessary calculation
            ntimes = 5;
            if k > ntimes
                d_rmse_n = diff(rmse_n(k-ntimes:k,1));
                if all(d_rmse_n>0)
                    % increased <ntimes> times in a row
                    % disp('breaking')
                    break
                end
            end

            % Find the value of smoothing parameter that minimizes
            % in-sample 1-step ahead prediction error
            minrmse = min(rmses);
            alpha_optimum = alphar(rmses == minrmse);

            % Produce some preliminary outputs
            out.train_minrmse = minrmse;
            out.train_alpha = alpha_optimum;

        end
    else
        % do a descent
        % fits a quadratic to available points

        % (1) use alpha = 0.01, 0.1, 0.5, 0.8
%         alphar = [0.1, 0.2, 0.5, 0.8, 0.9];
        alphar = linspace(0.1,0.9,5);
        rmses = zeros(4,1);

        for k = 1:length(alphar)
            a = alphar(k);

            xf = SUB_fit_exp_smooth(xtrain,a);

            % Issue forecasts
            fore = xf(3:end);
            orig = xtrain(3:end);
            rmses(k) = sqrt(mean((fore - orig).^2)); % compute rmse
        end

        % fit quadratic to set alpha
        [sort_rmses, ix] = sort(rmses);
        rkeep = ix(1:3); % fit on 3 points closest to minimum
        p = polyfit(alphar(rkeep)',rmses(rkeep),2);
        aar = (0:0.005:1);
        y = polyval(p,aar);
%         plot(aar,y,':k'); hold on;
%         plot(alphar,rmses,'or');
%         plot(alphar(rkeep),rmses(rkeep),'*m'); hold off
        alphamin = -p(2)/(2*p(1));
        out.alphamin_1 = alphamin;
        out.p1_1 = abs(p(1)); % concavity
        out.cup_1 = sign(p(1));

        if p(1) < 0 % concave down -- it's looking at a maximum
            % weird case
            if y(1) < y(end);
                alphamin = 0.01;
            else
                alphamin = 1;
            end
        else
            % Search again around this
            alphar = linspace(alphamin-0.1,alphamin+0.1,5);
            if any(alphar <= 0)
                alphar = linspace(0.01,max(alphamin,0)+0.1,5);
            elseif any(alphar >= 1)
                alphar = linspace(min(alphamin,1)-0.1,1,5);
            end

            for k = 1:length(alphar)
                a = alphar(k);

                xf = SUB_fit_exp_smooth(xtrain,a);

                % Issue forecasts
                fore = xf(3:end);
                orig = xtrain(3:end);
                rmses(k) = sqrt(mean((fore - orig).^2)); % compute rmse

            end

            % fit quadratic to set alpha
            p = polyfit(alphar',rmses,2);
%             aar = 0:0.005:1;
%             y = polyval(p,aar);
%             plot(aar,y,':k'); hold on;
%             plot(alphar,rmses,'or'); hold off

            if p(1) < 0
                alphamin = alphar(rmses == min(rmses));
                % This is quite bad -- the first step didn't find a local
                % minimum...
            else % minimum of quadratic fit
                alphamin = -p(2)/(2*p(1));
                if alphamin > 1, alphamin = 1; end
                if alphamin <= 0, alphamin = 0.01; end
            end

        end
    end

    out.alphamin = alphamin;
    alpha = alphamin;
end

if isnan(alpha)
    error('Alpha is a NaN?!')
end

%% (2) Fit to the whole time series

% % Plot in-sample error as a function of smoothing parameter
% figure(1);
% plot(alphar,rmses,'*');
% xlabel('\alpha - Smoothing parameter');
% ylabel('RMSE');

%Plot original time series and smoothed data using optimum values
y = SUB_fit_exp_smooth(x,alpha);

yp = y(3:N); % predicted
xp = x(3:N); % original
e = yp-xp; % residuals
% in_sample_error = sqrt(mean((yp-xp).^2));
% out.insamplermse = in_sample_error;


% Use MF_ResidualAnalysis on the residuals
% 1) Get statistics on residuals
residout = MF_ResidualAnalysis(e);

% convert these to local outputs in quick loop
fields = fieldnames(residout);
for k = 1:length(fields)
    out.(fields{k}) = residout.(fields{k});
end

% t=1:length(yp);

if doPlot
    figure('color','w'); box('on')
    plot(t,x(3:N),'b',t,y(3:N),'k');
    legend('Obs', 'Fit');
    xlabel('Time');
    ylabel('Amplitude');
end

% ------------------------------------------------------------------------------
function xf = SUB_fit_exp_smooth(x,a)
    % Iterate over rolling window
    ntrain = length(x);
    xf = zeros(ntrain,1);

    for ii = 2:ntrain-1
        s = zeros(ntrain,1);
        s(1) = mean(x(1:ii-1));

        % Loop to smooth data within the window size
        for jj = 2:ii
            s(jj) = a*x(jj) + (1-a)*s(jj-1);
        end

        % S(t) = Xf(t) is forecasted value for X(t+1)
        xf(ii+1,1) = s(ii);
    end
end


end
