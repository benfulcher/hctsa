% DV_dynsys_force
% 
% Couples the values of the time series to a given dynamical system. The input
% time series forces a particle in the given potential well.
% 
% The time series contributes to a forcing term on a simulated particle in a:
% (i) quartic double-well potential with potential energy V(x) = x^4/4 - alpha^2
% x^2/2, or a force F(x) = -x^3 + alpha^2 x, or
% (ii) a sinusoidal potential with V(x) = -cos(x/alpha), or F(x) = sin(x/alpha)/alpha.
%
% INPUTS:
% y, the input time series
% whatpot, the potential function to simulate:
%               'dblwell' (a double well potential function)
%               'sine' (a sinusoidal potential function)
% params, the parameters for simulation
%            should be in the form: params = [alpha, kappa, deltat]
% The double-well potential has three parameters:
% alpha controls the positions of the wells,
% kappa controls the coefficient of friction,
% and the step size of the simulation is Delta t.
% 
% The sinusoidal potential has parameters:
% alpha controls the period of oscillations in the potential,
% kappa is the coefficient of friction,
% Delta t is the time step.

% 
% Outputs are statistics summarizing the trajectory of the simulated particle,
% including its mean, the range, proportion positive, proportion of times it
% crosses zero, its autocorrelation, final position, and standard deviation.

function out = DV_dynsys_force(y,whatpot,params)
% Ben Fulcher, September 2009

% Check inputs
if nargin < 2 || isempty(whatpot)
    whatpot = 'dblwell'; % by default
end
if nargin < 3 || isempty(params)
    % default parameters
    switch whatpot
    case 'dblwell'
        params = [2, 0.1, 0.1];
    case 'sine'
        params = [1, 1, 1];
    otherwise
        error('Unknown system ''%s''', whatpot);
    end
end

doplot = 0; % plot results
N = length(y); % length of the time series

alpha = params(1);
kappa = params(2);
deltat = params(3); % time step

% Specify the potential function
switch whatpot
case 'sine'
    V = @(x) -cos(x/alpha);
    F = @(x) sin(x/alpha)/alpha;
case 'dblwell'
    F = @(x) -x.^3 + alpha^2*x; % the double well function (the force from a double well potential)
    V = @(x) x.^4/4 - alpha^2*x.^2/2;
otherwise
    error('Unknown potential function ''%s'' specified', whatpot);
end
    
x = zeros(N,1); % Position
v = zeros(N,1); % Velocity

for i = 2:N
    x(i) = x(i-1) + v(i-1)*deltat+(F(x(i-1))+y(i-1)-kappa*v(i-1))*deltat^2;
    v(i) = v(i-1) + (F(x(i-1))+y(i-1)-kappa*v(i-1))*deltat;
end

if doplot
    switch whatpot
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
    out = NaN; % not suitable for this time series
    return
end

%% Output some basic features of the trajectory
out.mean = mean(x); % mean
out.median = median(x); % median
out.std = std(x); % standard deviation
out.range = range(x); % range
out.proppos = sum(x > 0)/N; % proportion positive
out.pcross = sum((x(1:end-1)).*(x(2:end)) < 0)/(N-1); % n crosses middle
out.ac1 = abs(CO_autocorr(x,1)); % magnitude of autocorrelation at lag 1
out.ac10 = abs(CO_autocorr(x,10)); % magnitude of autocorrelation at lag 10
out.ac50 = abs(CO_autocorr(x,50)); % magnitude of autocorrelation at lag 50
out.tau = CO_fzcac(x); % first zero crossing of the autocorrelation function
out.finaldev = abs(x(end)); % final position

% A couple of additional outputs for double well:
if strcmp(whatpot,'dblwell')
    % number of times the trajectory crosses the middle of the upper well
    out.pcrossup = sum((x(1:end-1)-alpha).*(x(2:end)-alpha) < 0)/(N-1);
    % number of times the trajectory crosses the middle of the lower well
    out.pcrossdown = sum((x(1:end-1)+alpha).*(x(2:end)+alpha) < 0)/(N-1);
end


end