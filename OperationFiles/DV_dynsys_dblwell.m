function out=DV_dynsys_dblwell(y,params)
% Drives a dynamical system: double well potential using the vector of drives
% given in y
% Ben Fulcher September 2009

N = length(y);

alpha = params(1);
kappa = params(2);
deltat = params(3);

% V = @(x) x.^4/4-alpha^2*x.^2/2;
F = @(x) -x.^3+alpha^2*x;

x=zeros(N,1); % position
v=zeros(N,1); % velocity

for i=2:N
    x(i)=x(i-1)+v(i-1)*deltat+(F(x(i-1))+y(i-1)-kappa*v(i-1))*deltat^2;
    v(i)=v(i-1)+(F(x(i-1))+y(i-1)-kappa*v(i-1))*deltat;
end

% out=x;
% hold on;
% plot(-100:0.1:100,F(-100:0.1:100),'k')
% plot(x,V(x),'or')
% plot(x)
% hold off

%% OUTPUTS!
% features of the trajectory
if isnan(x(end)) || abs(x(end))>1E10
    out.mean=NaN; out.median=NaN; out.range=NaN; out.proppos=NaN; out.pcross=NaN;
    out.pcrossup=NaN; out.pcrossdown=NaN; out.ac1=NaN; out.ac10=NaN; out.ac50=NaN;
    out.tau=NaN; out.finaldev=NaN; out.std=NaN;
else
    out.mean = mean(x);
    out.median = median(x);
    out.range = range(x);
    out.proppos = length(find(x>0))/N; % proportion positive
    out.pcross = length(find((x(1:end-1)).*(x(2:end))<0))/(N-1); % crosses middle
    out.pcrossup = length(find((x(1:end-1)-alpha).*(x(2:end)-alpha)<0))/(N-1); % crosses upper well middle
    out.pcrossdown = length(find((x(1:end-1)+alpha).*(x(2:end)+alpha)<0))/(N-1); % crosses lower well middle
    out.ac1 = abs(CO_autocorr(x,1));
    out.ac10 = abs(CO_autocorr(x,10));
    out.ac50 = abs(CO_autocorr(x,50));
    out.tau = CO_fzcac(x);
    out.finaldev = abs(x(end));
    out.std = std(x);
end
% xth=SUB_coursegrain(x,2);



end