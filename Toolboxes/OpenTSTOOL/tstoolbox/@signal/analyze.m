function analyze(s, maxdim)

%tstoolbox/@signal/analyze
%   Syntax:
%     * analyze(s, maxdim)
%
%   Input arguments:
%     * maxdim - analyze will not use a dimension higher than this limit
%
%   Try to do a automatic analysis procedure of a time series. The time
%   series is embedded using the first zero of the auto mutual information
%   function for the delay time.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt



narginchk(2,2);

if (ndim(s) > 1) | (~isreal(data(s)))
    error('Input signal must be a scalar time series');
end

N = dlens(s,1); 		% length of time series
S = center(setaxis(s, 1, achse));	% set an x-axis without any unit, sampling rate etc.

tau1 = ceil(firstzero(acf(S)))
tau2 = firstmin(amutual(S, tau1 * 4))

disp(['Estimated delay time : ' num2str(tau2)]);

e = embed(s, maxdim, tau2);

E = data(e);

Nuse = N - tau1 * 50

view(signal(crossprediction(E(1:Nuse,:), E, 1:7:Nuse, tau1*20, ceil(N/1000), tau2)));
