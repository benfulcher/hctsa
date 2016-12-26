function out = PH_ForcePotential(y,whatPotential,params)
% PH_ForcePotential   Couples the values of the time series to a dynamical system
%
% The input time series forces a particle in the given potential well.
%
% The time series contributes to a forcing term on a simulated particle in a:
%
% (i) Quartic double-well potential with potential energy V(x) = x^4/4 - alpha^2
%           x^2/2, or a force F(x) = -x^3 + alpha^2 x
%
% (ii) Sinusoidal potential with V(x) = -cos(x/alpha), or F(x) = sin(x/alpha)/alpha
%
%---INPUTS:
% y, the input time series
%
% whatPotential, the potential function to simulate:
%               (i) 'dblwell' (a double well potential function)
%               (ii) 'sine' (a sinusoidal potential function)
%
% params, the parameters for simulation, should be in the form:
%                   params = [alpha, kappa, deltat]
%
%           (i) The double-well potential has three parameters:
%               * alpha controls the relative positions of the wells,
%               * kappa is the coefficient of friction,
%               * deltat sets the time step for the simulation.
%
%           (ii) The sinusoidal potential also has three parameters:
%               * alpha controls the period of oscillations in the potential
%               * kappa is the coefficient of friction,
%               * deltat sets the time step for the simulation.
%
%---OUTPUTS:
% Statistics summarizing the trajectory of the simulated particle,
% including its mean, the range, proportion positive, proportion of times it
% crosses zero, its autocorrelation, final position, and standard deviation.

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

doPlot = 0; % plot results
N = length(y); % length of the time series

% ------------------------------------------------------------------------------


alpha = params(1);
kappa = params(2);
deltat = params(3); % time step

% Specify the potential function
switch whatPotential
case 'sine'
    V = @(x) -cos(x/alpha);
    F = @(x) sin(x/alpha)/alpha;
case 'dblwell'
    F = @(x) -x.^3 + alpha^2*x; % the double well function (the force from a double well potential)
    V = @(x) x.^4/4 - alpha^2*x.^2/2;
otherwise
    error('Unknown potential function ''%s'' specified', whatPotential);
end

x = zeros(N,1); % Position
v = zeros(N,1); % Velocity

for i = 2:N
    x(i) = x(i-1) + v(i-1)*deltat+(F(x(i-1))+y(i-1)-kappa*v(i-1))*deltat^2;
    v(i) = v(i-1) + (F(x(i-1))+y(i-1)-kappa*v(i-1))*deltat;
end

if doPlot
    switch whatPotential
    case 'dblwell'
        figure('color','w'); hold on;
        plot(-100:0.1:100, F(-100:0.1:100), 'k') % plot the potential
        plot(x,V(x),'or')
        plot(x)

    case 'sine'
        figure('color','w');
        subplot(3,1,1); plot(y,'k'); title('Time series -> drive')
        subplot(3,1,2); plot(x,'k'); title('Simulated particle position')
        subplot(3,1,3); box('on'); hold on;
        plot(min(x):0.1:max(x),V(min(x):0.1:max(x)),'k')
        plot(x,V(x),'.r')
    end
end

% Check trajectory didn't blow out:
if isnan(x(end)) || abs(x(end)) > 1E10
    fprintf(1,'Trajectory blew out!\n');
    out = NaN; return % not suitable for this time series
end

% ------------------------------------------------------------------------------
%% Output some basic features of the trajectory
% ------------------------------------------------------------------------------
out.mean = mean(x); % mean
out.median = median(x); % median
out.std = std(x); % standard deviation
out.range = range(x); % range
out.proppos = sum(x > 0)/N; % proportion positive
out.pcross = sum((x(1:end-1)).*(x(2:end)) < 0)/(N-1); % n crosses middle
out.ac1 = abs(CO_AutoCorr(x,1,'Fourier')); % magnitude of autocorrelation at lag 1
out.ac10 = abs(CO_AutoCorr(x,10,'Fourier')); % magnitude of autocorrelation at lag 10
out.ac50 = abs(CO_AutoCorr(x,50,'Fourier')); % magnitude of autocorrelation at lag 50
out.tau = CO_FirstZero(x,'ac'); % first zero crossing of the autocorrelation function
out.finaldev = abs(x(end)); % final position

% A couple of additional outputs for double well:
if strcmp(whatPotential,'dblwell')
    % number of times the trajectory crosses the middle of the upper well
    out.pcrossup = sum((x(1:end-1)-alpha).*(x(2:end)-alpha) < 0)/(N-1);
    % number of times the trajectory crosses the middle of the lower well
    out.pcrossdown = sum((x(1:end-1)+alpha).*(x(2:end)+alpha) < 0)/(N-1);
end


end
