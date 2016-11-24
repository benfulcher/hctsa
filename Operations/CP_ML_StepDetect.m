function out = CP_ML_StepDetect(y,method,params)
% CP_ML_StepDetect      Analysis of discrete steps in a time series.
%
% Gives information about discrete steps in the signal, using the function
% l1pwc from Max A. Little's step detection toolkit.
%
% cf.,
% "Sparse Bayesian Step-Filtering for High-Throughput Analysis of Molecular
% Machine Dynamics", Max A. Little, and Nick S. Jones, Proc. ICASSP (2010)
%
% "Steps and bumps: precision extraction of discrete states of molecular machines"
% M. A. Little, B. C. Steel, F. Bai, Y. Sowa, T. Bilyard, D. M. Mueller,
% R. M. Berry, N. S. Jones. Biophysical Journal, 101(2):477-485 (2011)
%
% Software available at: http://www.maxlittle.net/software/index.php
%
%---INPUTS:
% y, the input time series
%
% method, the step-detection method:
%           (i) 'kv': Kalafut-Visscher
%                 cf. The algorithm described in:
%                 Kalafut, Visscher, "An objective, model-independent method for
%                 detection of non-uniform steps in noisy signals", Comp. Phys.
%                 Comm., 179(2008), 716-723.
%
%           (ii) 'l1pwc': L1 method
%                 This code is based on code originally written by Kim et al.:
%                 "l_1 Trend Filtering", S.-J. Kim et al., SIAM Review 51, 339
%                 (2009).
%
% params, the parameters for the given method used:
%           (i) 'kv': (no parameters required)
%           (ii) 'l1pwc': params = lambda
%
%---OUTPUTS:
% Statistics on the output of the step-detection method, including the intervals
% between change points, the proportion of constant segments, the reduction in
% variance from removing the piece-wise constants, and stationarity in the
% occurrence of change points.

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
% Preliminaries, check inputs:
% ------------------------------------------------------------------------------
N = length(y); % time-series length

if nargin < 2 || isempty(method)
    fprintf(1,'Using Kalafut-Visscher step detection by default\n');
    method = 'kv';
end

% ------------------------------------------------------------------------------
% Do the step detection:
% ------------------------------------------------------------------------------

switch method
    case 'kv'
        % ------------------------------------------------------------------------------
        %% Kalafut-Visscher
        % ------------------------------------------------------------------------------

        % (1) Do the step detection
        [steppedy, steps] = ML_kvsteps(y);

        % Put in chpts form: a vector specifying indicies of starts of
        % constant runs.
        if length(steps) == 2
            chpts = 1;
        else
            chpts = [1; steps(2:end-1)+1];
        end

    % case 'ck'
    %     % ------------------------------------------------------------------------------
    %     %% Chung-Kennedy
    %     % ------------------------------------------------------------------------------
    %     % The algorithm is described in:
    %     % S.H. Chung, R.A. Kennedy (1991), "Forward-backward non-linear filtering
    %     % technique for extracting small biological signals from noise",
    %     % J. Neurosci. Methods. 40(1):71-86.
    %     % It is quite slow...
    %     % And not supported!
    %
    %     % Inputs:
    %     %  y - Input signal
    %     %  K - Maximum forward/backward moving average filter length (samples)
    %     %  M - Prediction error analysis window size (samples)
    %     %  p - Positive scaling of prediction error
    %     % Outputs:
    %     %  x - Step-filtered output signal
    %
    %     % Set defaults, params should be [K,M,p]
    %     if nargin < 3
    %         params = [];
    %     end
    %     if length(params) >= 1
    %         K = params(1);
    %     else
    %         K = 1/20; % 1/20th of time series length
    %     end
    %     if K < 1
    %         K = floor(N*K);
    %     end
    %     if length(params) >= 2
    %         M = params(2);
    %     else
    %         M = 1/10; % 1/10th the time series length
    %     end
    %     if M < 1
    %         M = floor(N*M);
    %     end
    %     if length(params) >= 3
    %         p = params(3);
    %     else
    %         p = 10;
    %     end
    %     steppedy = ML_ckfilter(y, K, M, p);

    case 'l1pwc'
        % ------------------------------------------------------------------------------
        % Based around code originally written by
        % S.J. Kim, K. Koh, S. Boyd and D. Gorinevsky. If you use this code for
        % your research, please cite:
        % M.A. Little, Nick S. Jones (2010)
        % "Sparse Bayesian Step-Filtering for High-Throughput Analysis of Molecular
        % Machine Dynamics", in 2010 IEEE International Conference on Acoustics,
        % Speech and Signal Processing, 2010, ICASSP 2010 Proceedings.
        % ------------------------------------------------------------------------------

        % Input arguments:
        % - y          Original signal to denoise, size N x 1.
        % - lambda     A vector of positive regularization parameters, size L x 1.
        %              TVD will be applied to each value in the vector.
        % - display    (Optional) Set to 0 to turn off progress display, 1 to turn
        %              on. If not specifed, defaults to progress display on.
        % - stoptol    (Optional) Precision as determined by duality gap tolerance,
        %              if not specified, defaults to 1e-3.
        % - maxiter    (Optional) Maximum interior-point iterations, if not
        %              specified defaults to 60.
        %
        % Output arguments:
        % - x          Denoised output signal for each value of lambda, size N x L.
        % - E          Objective functional at minimum for each lambda, size L x 1.
        % - s          Optimization result, 1 = solved, 0 = maximum iterations
        %              exceeded before reaching duality gap tolerance, size L x 1.
        % - lambdamax  Maximum value of lambda for the given y. If
        %              lambda >= lambdamax, the output is the trivial constant
        %              solution x = mean(y).

        % Set defaults, params should be [lambda]
        if nargin < 3
            params = [];
        end
        if length(params) >= 1
            lambda = params(1);
        else
            lambda = 10; % higher lambda --> less steps
        end
        if lambda < 1 % specify as a proportion of lambdamax
            lambda = ML_l1pwclmax(y)*lambda;
        end

        % Run the code
        [steppedy, E, s, lambdamax] = ML_l1pwc(y, lambda, 0); % use defaults for stoptol and maxiter

        % round to this precision to remove numberical flucuations of order
        % less than 1e-4
        steppedy = round(steppedy*1e4)/1e4;

        % Return outputs
        out.E = E;
        out.s = s; % for some parameter values, this is 1
        out.lambdamax = lambdamax;

        % Get step indicies from steppedy
        % these give the index of the start of each run
        whch = find(diff(steppedy)~=0);
        if ~isempty(whch)
            chpts = [1;whch+1];
        else
            chpts = 1; % no changes
        end

    otherwise
        error('Unknown step detection method ''%s''',method);
end

% ------------------------------------------------------------------------------
% Outputs common to all step detection methods:
% ------------------------------------------------------------------------------
% requires: chpts -- a vector of indicies for changes in the time series
%           steppedy -- an (Nx1) vector specifying the new stepped time
%                       series

% Intervals -- of change
chints = diff([chpts; N]);

% Number of constant segments per sample
out.nsegments = length(chpts)/N; % will be 1 if there are no changes
% How much reduces variance
out.rmsoff = std(y) - std(y-steppedy);
% Reduces variance per step
out.rmsoffpstep = (out.rmsoff)/(length(chpts));
% Ratio of number of steps in first half of time series to second half
sum1 = (sum(chpts < N/2)-1);
sum2 = sum(chpts >= N/2);
if (sum2 > 0) && (sum1 > 0)
    if sum2 > sum1
        out.ratn12 = sum1/sum2;
    else
        out.ratn12 = sum2/sum1;
    end
else
    out.ratn12 = 0;
end
out.diffn12 = abs(sum1-sum2)/length(chpts);

% Proportion of really short steps:
out.pshort_3 = sum(chints <= 3)/N;
% Mean interval between steps:
out.meanstepint = mean(chints)/N;
% Mean interval greater than 3 samples, per sample:
out.meanstepintgt3 = mean(chints(chints>3))/N;
% Mean error on step interval distribution:
out.meanerrstepint = std(chints)/sqrt(length(chints));
% Maximum step interval:
out.maxstepint = max(chints)/N;
% Minimum step interval:
out.minstepint = min(chints)/N;
% Median step interval:
out.medianstepint = median(chints)/N;


end
