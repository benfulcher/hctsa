function out = CO_basicrecurf(y,tau)
% Analyzes the recurrence plot of y(i) against y(i-tau)
% INPUTS:
% y is the z-scored time series
% tau is the time lag (can be 'tau' to set the time lag the first zero crossing of the autocorrelation function)
% Returns basic measures of structure in this two-dimensional recurrence space
% Ben Fulcher September 2009

doplot = 0;

if strcmp(tau,'tau')
	% Make tau the first zero crossing of the autocorrelation function
    tau = CO_fzcac(y);
end

xt = y(1:end-tau); % part of the time series
xtp = y(1+tau:end); % time-lagged time series
N = length(y) - tau; % Length of each time series subsegment


% Points in a thick bottom-left -- top-right diagonal
out.updiag01 = sum(abs(xtp-xt) < 0.1)/N;
out.updiag05 = sum(abs(xtp-xt) < 0.5)/N;

% Points in a thick bottom-right -- top-left diagonal
out.downdiag01 = sum(abs(xtp+xt) < 0.1)/N;
out.downdiag05 = sum(abs(xtp+xt) < 0.5)/N;

% Ratio of these
out.ratdiag01 = out.updiag01/out.downdiag01;
out.ratdiag05 = out.updiag05/out.downdiag05;

% In a thick parabola concave up
out.parabup01 = sum(abs(xtp-xt.^2) < 0.1)/N;
out.parabup05 = sum(abs(xtp-xt.^2) < 0.5)/N;

% In a thick parabola concave down
out.parabdown01 = sum(abs(xtp+xt.^2) < 0.1)/N;
out.parabdown05 = sum(abs(xtp+xt.^2) < 0.5)/N;

% In a thick parabola concave up, shifted up 1
out.parabup01_1 = sum(abs(xtp-(xt.^2+1)) < 0.1)/N;
out.parabup05_1 = sum(abs(xtp-(xt.^2+1)) < 0.5)/N;

% In a thick parabola concave down, shifted up 1
out.parabdown01_1 = sum(abs(xtp+(xt.^2-1)) < 0.1)/N;
out.parabdown05_1 = sum(abs(xtp+(xt.^2-1)) < 0.5)/N;

% In a thick parabola concave up, shifted down 1
out.parabup01_n1 = sum(abs(xtp-(xt.^2-1)) < 0.1)/N;
out.parabup05_n1 = sum(abs(xtp-(xt.^2-1)) < 0.5)/N;

% In a thick parabola concave down, shifted down 1
out.parabdown01_n1 = sum(abs(xtp+(xt.^2+1)) < 0.1)/N;
out.parabdown05_n1 = sum(abs(xtp+(xt.^2+1)) < 0.5)/N;

% RINGS (points within a radius range)
out.ring1_01 = sum(abs(xtp.^2+xt.^2-1) < 0.1)/N;
out.ring1_02 = sum(abs(xtp.^2+xt.^2-1) < 0.2)/N;
out.ring1_05 = sum(abs(xtp.^2+xt.^2-1) < 0.5)/N;

% CIRCLES (points inside a given circular boundary)
out.incircle_01 = sum(xtp.^2+xt.^2 < 0.1)/N;
out.incircle_02 = sum(xtp.^2+xt.^2 < 0.2)/N;
out.incircle_05 = sum(xtp.^2+xt.^2 < 0.5)/N;
out.incircle_1 = sum(xtp.^2+xt.^2 < 1)/N;
out.incircle_2 = sum(xtp.^2+xt.^2 < 2)/N;
out.incircle_3 = sum(xtp.^2+xt.^2 < 3)/N;
out.medianincircle = median([out.incircle_01, out.incircle_02, out.incircle_05 ...
                            out.incircle_1, out.incircle_2, out.incircle_3]);
out.stdincircle = std([out.incircle_01, out.incircle_02, out.incircle_05 ...
                        out.incircle_1, out.incircle_2, out.incircle_3]);


if doplot
    figure('color','w');
    plot(xt,xtp,'.k');
    hold on
    r = (xtp.^2+xt.^2 < 0.2);
    plot(xt(r),xtp(r),'.g')
end

end