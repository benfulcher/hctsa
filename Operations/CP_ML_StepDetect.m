function out = CP_ML_StepDetect(y, method, params)
% CP_ML_StepDetect   Fit a staircase of constant levels to the series and describe its steps.
%
% Finds discrete steps in the series using the l1pwc function from Max A.
% Little's step-detection toolkit. With method 'l1pwc', the series is replaced by
% the piecewise-constant signal s that minimizes
%       (1/2) sum_t (y_t - s_t)^2 + lambda * sum_t |s_{t+1} - s_t|
% (total-variation denoising), where lambda is a penalty on the size of each step:
% a larger lambda gives fewer, larger steps. The outputs describe the fit (its
% objective value, the variance it removes) and the constant segments it finds
% (how many per sample, how long, how variable in length, and how evenly the
% change points are spread across the first and second halves of the series).
%
% ---INPUTS:
% y, the input time series
%
% method, the step-detection method:
%       (i) 'kv': Kalafut-Visscher (not used by hctsa),
%       (ii) 'l1pwc': L1 method (total-variation denoising). Based on code by
%            Kim et al. for l_1 trend filtering; here as implemented by Max A.
%            Little. The default if no method is given is 'kv'.
%
% params, the parameters for the method:
%       (i) 'kv': (no parameters required)
%       (ii) 'l1pwc': params = lambda, the penalty on each step's size
%            (default 10). For a z-scored series a fixed lambda gives a
%            segmentation that does not depend on the series length. A value
%            lambda < 1 is instead taken as a proportion of lambdamax (the
%            smallest lambda that gives no steps), which grows with the length,
%            as ~sqrt(N) for a short-memory process.
%
% ---OUTPUTS:
% For both methods:
% nsegments, the number of constant segments per sample (1/N if there are no steps)
% rmsoff, the reduction in standard deviation: std(y) - std(y - fit)
% rmsoffpstep, rmsoff divided by the number of constant segments
% ratn12, ratio of the number of change points in the first half of the series to
%       the number in the second half (smaller over larger; 0 if either is zero)
% diffn12, absolute difference between the number of change points in the two
%       halves, as a proportion of the number of segments
% pshort_3, number of segments of 3 samples or fewer, per sample
% meanstepintgt3, mean length (in samples) of the segments longer than 3 samples
% cvstepint, coefficient of variation of the segment lengths
% medianstepint, median segment length (in samples)
% For method 'l1pwc' only:
% E, the value of the objective above at its minimum, per sample (E/N)
% s, whether the solver converged (1) or hit its maximum iterations (0)
% lambdamax, the smallest lambda that gives no steps, divided by sqrt(N)
%
% ---REFERENCES:
% Max A. Little and Nick S. Jones, "Sparse Bayesian Step-Filtering for
% High-Throughput Analysis of Molecular Machine Dynamics", Proc. ICASSP (2010).
%
% M. A. Little, B. C. Steel, F. Bai, Y. Sowa, T. Bilyard, D. M. Mueller,
% R. M. Berry, N. S. Jones, "Steps and bumps: precision extraction of discrete
% states of molecular machines", Biophysical Journal 101(2): 477-485 (2011).
%
% Kalafut and Visscher, "An objective, model-independent method for detection of
% non-uniform steps in noisy signals", Comp. Phys. Comm. 179, 716-723 (2008).
%
% S.-J. Kim et al., "l_1 Trend Filtering", SIAM Review 51, 339 (2009).
%
% Software available at: http://www.maxlittle.net/software/index.php

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

% -------------------------------------------------------------------------------
% Check inputs:
if nargin < 2 || isempty(method)
	fprintf(1, 'Using Kalafut-Visscher step detection by default\n');
	method = 'kv';
end

% -------------------------------------------------------------------------------
% Preliminaries:
doPlot = false; % whether to plot outputs
N = length(y); % time-series length

% ------------------------------------------------------------------------------
% Do the step detection:
% ------------------------------------------------------------------------------
switch method
	case 'kv'
		% Kalafut-Visscher step detection
		[steppedy, steps] = ML_kvsteps(y);

		% Put in chpts form: a vector specifying indicies of starts of
		% constant runs.
		if length(steps) == 2
			chpts = 1;
		else
			chpts = [1; steps(2:end - 1) + 1];
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
		% Outputs:
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
			lambda = ML_l1pwclmax(y) * lambda;
		end

		% Run the code
		[steppedy, E, s, lambdaMax] = ML_l1pwc(y, lambda, 0); % use defaults for stoptol and maxiter

		% Round to remove numberical flucuations of order less than 1e-4
		steppedy = round(steppedy * 1e4) / 1e4;

		% Compute outputs specific to this method:
		out.E = E / N; % energy per sample (E sums over the series, so it scales with length)
		out.s = s; % for some parameter values, this is 1
		out.lambdamax = lambdaMax / sqrt(N); % (lambdamax itself grows as ~sqrt(N))

		% Get step indicies from steppedy
		% these give the index of the start of each run
		whch = find(diff(steppedy) ~= 0);
		if ~isempty(whch)
			chpts = [1; whch + 1];
		else
			chpts = 1; % no changes
		end
	otherwise
		error('Unknown step detection method ''%s''', method);
end

% -------------------------------------------------------------------------------
% Plot computed change points onto the time series:
if doPlot
	f = figure('color', 'w');
	hold('on')
	plot(y, 'k')
	plot(chpts, y(chpts), 'or')
end

% -------------------------------------------------------------------------------
% Outputs common to all step detection methods:
% requires: chpts -- a vector of indicies for changes in the time series
%           steppedy -- an (Nx1) vector specifying the new stepped time
%                       series

numChangePoints = length(chpts);

% Intervals -- of change (the length of each constant segment)
chints = diff([chpts; N + 1]);

% Number of constant segments per sample
out.nsegments = numChangePoints / N; % will be 1 if there are no changes

% How much reduces variance
out.rmsoff = std(y) - std(y - steppedy);

% Reduces variance per segment (numChangePoints counts the segment starting at
% sample 1, so it is at least 1 and this is defined even when there are no steps)
out.rmsoffpstep = out.rmsoff / numChangePoints;

% Ratio of number of steps in first half of time series to second half
sum1 = sum(chpts < N / 2) - 1; % (exclude the chpt that's always sitting at 1)
sum2 = sum(chpts >= N / 2);
if (sum2 > 0) && (sum1 > 0)
	if sum2 > sum1
		out.ratn12 = sum1 / sum2;
	else
		out.ratn12 = sum2 / sum1;
	end
else
	out.ratn12 = 0;
end
% Difference between number of change points in first and second half of data
% as a proportion of total number of change points
% Maximal (1) when all change points are in one half of data
% Minimal (0) when same number of change points in both halves
out.diffn12 = abs(sum1 - sum2) / numChangePoints;

% Proportion of really short steps:
out.pshort_3 = sum(chints <= 3) / N;
% Step intervals, in samples:
% (With a fixed lambda these do not grow with the series length. The mean
% interval is 1/nsegments, and so is not given separately.)
% Mean interval greater than 3 samples:
out.meanstepintgt3 = mean(chints(chints > 3));
% Coefficient of variation of the step intervals:
out.cvstepint = std(chints) / mean(chints);
% Median step interval:
out.medianstepint = median(chints);

end
