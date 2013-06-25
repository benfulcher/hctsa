function out = DV_dynsys_dblwell(y,params)
% Drives a dynamical system: double well potential using the vector of drives
% The input time series is given as y as forces a particle in the double well potential
% Parameters are in the form params = [alpha, kappa, deltat]
% Ben Fulcher September 2009

doplot = 0; % plot results

N = length(y);

alpha = params(1);
kappa = params(2);
deltat = params(3); % time step

F = @(x) -x.^3 + alpha^2*x; % the double well function (the force from a double well potential)

x = zeros(N,1); % Position
v = zeros(N,1); % Velocity

for i = 2:N
    x(i) = x(i-1) + v(i-1)*deltat+(F(x(i-1))+y(i-1)-kappa*v(i-1))*deltat^2;
    v(i) = v(i-1) + (F(x(i-1))+y(i-1)-kappa*v(i-1))*deltat;
end

if doplot
    figure('color','w'); hold on;
    plot(-100:0.1:100, F(-100:0.1:100), 'k') % plot the potential
    V = @(x) x.^4/4 - alpha^2*x.^2/2;
    plot(x,V(x),'or')
    plot(x)
    hold off
end

%% Output features of the trajectory
if isnan(x(end)) || abs(x(end)) > 1E10
    error('Error simulating time series forcing a double-well potential')
    % out = NaN;
    % out.mean = NaN; out.median = NaN; out.range = NaN; out.proppos = NaN; out.pcross = NaN;
    % out.pcrossup = NaN; out.pcrossdown = NaN; out.ac1 = NaN; out.ac10 = NaN; out.ac50 = NaN;
    % out.tau = NaN; out.finaldev = NaN; out.std = NaN;
else
    out.mean = mean(x);
    out.median = median(x);
    out.range = range(x);
    out.proppos = length(find(x>0))/N; % proportion positive
    out.pcross = length(find((x(1:end-1)).*(x(2:end))<0))/(N-1);
    % crosses middle
    out.pcrossup = length(find((x(1:end-1)-alpha).*(x(2:end)-alpha)<0))/(N-1); % crosses upper well middle
    out.pcrossdown = sum((x(1:end-1)+alpha).*(x(2:end)+alpha)<0)/(N-1);
    % crosses lower well middle
    
    out.ac1 = abs(CO_autocorr(x,1)); % magnitude of autocorrelation at lag 1
    out.ac10 = abs(CO_autocorr(x,10)); % magnitude of autocorrelation at lag 10
    out.ac50 = abs(CO_autocorr(x,50)); % magnitude of autocorrelation at lag 50
    out.tau = CO_fzcac(x); % first zero crossing of the autocorrelation function
    out.finaldev = abs(x(end));
    out.std = std(x);
end

% xth = SUB_coursegrain(x,2);

end