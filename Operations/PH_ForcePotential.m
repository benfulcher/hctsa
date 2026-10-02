function out = PH_ForcePotential(y, whatPotential, params)
% PH_ForcePotential   Statistics of a simulated particle in a potential well, pushed by the time series.
%
% The time series is used as a time-varying force on a simulated particle (with
% position x and velocity v) that also feels a potential and friction. The
% particle is simulated for as many time steps as the series has points, with
% the force from the potential and the next value of the series added at each step:
%     x(t) = x(t-1) + v(t-1) dt + (F(x(t-1)) + y(t-1) - kappa v(t-1)) dt^2
%     v(t) = v(t-1) + (F(x(t-1)) + y(t-1) - kappa v(t-1)) dt
% The outputs are statistics of the trajectory x(t). The potentials are:
%
% (i) A quartic double-well potential with V(x) = x^4/4 - alpha^2 x^2/2, and so
%     force F(x) = -x^3 + alpha^2 x, with wells at x = +alpha and x = -alpha.
%
% (ii) A sinusoidal potential with V(x) = -cos(x/alpha), and so force
%     F(x) = -sin(x/alpha)/alpha.
%
% ---INPUTS:
% y, the input time series
% whatPotential, the potential function to simulate:
%       'dblwell': a double-well potential
%       'sine': a sinusoidal potential
%       Default: 'dblwell'.
% params, the parameters of the simulation, [alpha, kappa, deltat]:
%       alpha, for the double well, the position of the wells (+/-alpha); for the
%           sinusoid, sets the period of the oscillations in the potential
%       kappa, the coefficient of friction
%       deltat, the time step of the simulation
%       Defaults: [2, 0.1, 0.1] for 'dblwell' and [1, 1, 1] for 'sine'.
%
% ---OUTPUTS: statistics of the trajectory of the particle (a scalar NaN if the
% trajectory blows up, or ends beyond 1e10):
% mean, median, std, range, the mean, median, standard deviation and range of x
% proppos, the proportion of time steps at which x is positive
% pcross, the proportion of time steps at which x crosses zero
% pcrossup, pcrossdown (double well only), the proportions of time steps at which x
%       crosses the center of the upper (x = alpha) and lower (x = -alpha) well
% ac1, ac10, ac50, the magnitude (absolute value) of the autocorrelation of x at
%       lags 1, 10 and 50
% tau, the first zero-crossing of the autocorrelation function of x
% finaldev, the magnitude of the final position, |x(end)|
%
% ---NOTES:
% The update is the semi-implicit (symplectic) Euler scheme, x(t) = x(t-1) + v(t) dt
% with the new velocity v(t), which is the form written above. It is not the
% constant-acceleration formula (which has dt^2/2 and is not symplectic): that
% variant makes the double-well simulation blow up for many series and parameter sets.

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

% ------------------------------------------------------------------------------
% Check inputs and set defaults
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(whatPotential)
	whatPotential = 'dblwell'; % by default
end

if nargin < 3 || isempty(params)
	% default parameters
	switch whatPotential
		case 'dblwell'
			params = [2, 0.1, 0.1];
		case 'sine'
			params = [1, 1, 1];
		otherwise
			error('Unknown system ''%s''', whatPotential);
	end
end

doPlot = false; % plot results
N = length(y); % length of the time series

% ------------------------------------------------------------------------------

alpha = params(1);
kappa = params(2);
deltat = params(3); % time step

% Specify the potential function
switch whatPotential
	case 'sine'
		V = @(x) -cos(x / alpha);
		F = @(x) -sin(x / alpha) / alpha; % F = -dV/dx
	case 'dblwell'
		F = @(x) -x.^3 + alpha^2 * x; % the double well function (the force from a double well potential)
		V = @(x) x.^4 / 4 - alpha^2 * x.^2 / 2;
	otherwise
		error('Unknown potential function ''%s'' specified', whatPotential);
end

x = zeros(N, 1); % Position
v = zeros(N, 1); % Velocity

for i = 2:N
	x(i) = x(i - 1) + v(i - 1) * deltat + (F(x(i - 1)) + y(i - 1) - kappa * v(i - 1)) * deltat^2;
	v(i) = v(i - 1) + (F(x(i - 1)) + y(i - 1) - kappa * v(i - 1)) * deltat;
end

if doPlot
	switch whatPotential
		case 'dblwell'
			figure('color', 'w'); hold on;
			plot(-100:0.1:100, F(-100:0.1:100), 'k') % plot the potential
			plot(x, V(x), 'or')
			plot(x)

		case 'sine'
			figure('color', 'w');
			subplot(3, 1, 1); plot(y, 'k'); title('Time series -> drive')
			subplot(3, 1, 2); plot(x, 'k'); title('Simulated particle position')
			subplot(3, 1, 3); box('on'); hold on;
			plot(min(x):0.1:max(x), V(min(x):0.1:max(x)), 'k')
			plot(x, V(x), '.r')
	end
end

% Check trajectory didn't blow out:
if isnan(x(end)) || abs(x(end)) > 1E10
	fprintf(1, 'Trajectory blew out!\n');
	out = NaN;
	return % not suitable for this time series
end

% ------------------------------------------------------------------------------
%% Output some basic features of the trajectory
% ------------------------------------------------------------------------------
out.mean = mean(x); % mean
out.median = median(x); % median
out.std = std(x); % standard deviation
out.range = range(x); % range
out.proppos = sum(x > 0) / N; % proportion positive
out.pcross = sum((x(1:end - 1)) .* (x(2:end)) < 0) / (N - 1); % n crosses middle
out.ac1 = abs(CO_AutoCorr(x, 1, 'Fourier')); % magnitude of autocorrelation at lag 1
out.ac10 = abs(CO_AutoCorr(x, 10, 'Fourier')); % magnitude of autocorrelation at lag 10
out.ac50 = abs(CO_AutoCorr(x, 50, 'Fourier')); % magnitude of autocorrelation at lag 50
out.tau = CO_FirstCrossing(x, 'ac', 0, 'continuous'); % first zero-crossing of the autocorrelation function
out.finaldev = abs(x(end)); % final position

% A couple of additional outputs for double well:
if strcmp(whatPotential, 'dblwell')
	% number of times the trajectory crosses the middle of the upper well
	out.pcrossup = sum((x(1:end - 1) - alpha) .* (x(2:end) - alpha) < 0) / (N - 1);
	% number of times the trajectory crosses the middle of the lower well
	out.pcrossdown = sum((x(1:end - 1) + alpha) .* (x(2:end) + alpha) < 0) / (N - 1);
end

end
