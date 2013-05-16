% Simulate bacterial flagellar motor time series.
%
% Usage:
% [t, y, m] = bfmsim(dt, simdur, kappa, xi, shape, scale, ndwells, aperiod)
%
% Inputs
%    dt      - Seconds per sample
%    simdur  - Simulation duration in seconds
%    kappa   - Connecting spring constant
%    xi      - Stokes coefficient of friction
%    shape   - Dwell time gamma distribution shape parameter
%    scale   - Dwell time gamma distribution scale parameter
%    ndwells - Number of dwell locations in each rev
%    aperiod - Aperiodicity fraction of dwell locations, 0 = periodic
%
% Outputs
%    t  - Vector of simulation time points
%    y  - Observed motor time series
%    m  - Underlying motor position time series
%    mu - Dwell locations
%
% (c) Max Little, 2010. If you use this code for your research, please
% cite:
% "Steps and bumps: precision extraction of discrete states of molecular
% machines using physically-based, high-throughput time series analysis"
% Max A. Little et al., 2010, arXiv:1004.1234v1 [q-bio.QM]

function [t, y, m, mu] = bfmsim(dt, simdur, kappa, xi, shape, scale, ndwells, aperiod)

error(nargchk(8,8,nargin));

% Synthesise dwell locations
dtheta = 2*pi/ndwells;
mu = (0:(ndwells-1))*2*pi/ndwells+dtheta/2;
mu = mu + aperiod*dtheta*randn(1,ndwells);

% Calculate AR coefficient and noise variance
a = 1-kappa*dt/xi;
if (a > 0)
    % Here the relaxation slower than sample rate
    stdeps = sqrt(2*dt/xi);
else
    % Relaxation faster than sample rate - noise becomes uncorrelated
    stdeps = sqrt(kappa)*xi;
end

t = 0:dt:(simdur-dt);
N = length(t);
n = 1;
m = [];
while (length(m) < N)
    M = floor(gamrnd(shape, scale));
    m = [m; mu(n)*ones(M,1)];
    n = n + 1;
    if (n > ndwells)
        n = 1;
    end
end
m = m(1:N);
u = unwrap(m);
e = stdeps*randn(N,1);

if (a > 0)
    x = (1-a)*u + e;
    y = filter(1,[1 -a],x);
else
    y = u + e;
end
y = mod(y, 2*pi);
