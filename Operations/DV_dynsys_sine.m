function out = DV_dynsys_sine(y,params)
% Drives a dynamical system: double well potential using the vector of drives
% given in y
% Ben Fulcher September 2009

N = length(y);

alpha = params(1);
kappa = params(2);
deltat = params(3);

V = @(x) -cos(x/alpha);
F = @(x) sin(x/alpha)/alpha;

x = zeros(N,1); % position
v = zeros(N,1); % velocity

for i = 2:N
    x(i) = x(i-1)+v(i-1)*deltat+(F(x(i-1))+y(i-1)-kappa*v(i-1))*deltat^2;
    v(i) = v(i-1)+(F(x(i-1))+y(i-1)-kappa*v(i-1))*deltat;
end

% PLOT:
% figure('color','w');
% subplot(3,1,1); plot(y,'k'); title('time series -> drive')
% subplot(3,1,2); plot(x,'k'); title('simulated particle position')
% subplot(3,1,3); box('on'); hold on;
% plot(min(x):0.1:max(x),V(min(x):0.1:max(x)),'k')
% plot(x,V(x),'.r')

%% OUTPUTS!
% features of the trajectory
if isnan(x(end)) || abs(x(end)) > 1E10
    % trajectory blew out
    disp('DV_dynsys_sine: Trajectory blew out!');
    out = NaN;
%     out.mean=NaN; out.median=NaN; out.range=NaN; out.proppos=NaN; out.pcross=NaN;
%     out.ac1=NaN; out.ac10=NaN; out.ac50=NaN; out.tau=NaN; out.finaldev=NaN; out.std=NaN;
else
    out.mean = mean(x); % mean
    out.median = median(x); % median
    out.range = range(x); % range
    out.proppos = length(find(x>0))/N; % proportion positive
    out.pcross = length(find((x(1:end-1)).*(x(2:end))<0))/(N-1); % crosses middle
    out.ac1 = abs(CO_autocorr(x,1));
    out.ac10 = abs(CO_autocorr(x,10));
    out.ac50 = abs(CO_autocorr(x,50));
    out.tau = CO_fzcac(x);
    out.finaldev = abs(x(end));
    out.std = std(x);
end

% xth=SUB_coursegrain(x,2);



end