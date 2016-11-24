function out = SY_RangeEvolve(y)
% SY_RangeEvolve    How the time-series range changes across time.
%
% Measures of the range of the time series as a function of time,
% i.e., range(x_{1:i}) for i = 1, 2, ..., N, where N is the length of the time
% series.
%
%---INPUT:
% y, the time series
%
%---OUTPUTS: based on the dynamics of how new extreme events occur with time.

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

doPlot = 0; % plot outputs
N = length(y); % length of the time series
cums = zeros(N,1);

for i = 1:N
    cums(i) = range(y(1:i));
end
% cums=cums/range(y);

if doPlot
    figure('color','w');
    plot(cums);
end

fullr = range(y);

% Define an anonymous function for the number of unique entries in a vector, x:
lunique = @(x) length(unique(x));

out.totnuq = lunique(cums);

% How many of the unique extrema are in first <proportion> of time series?
cumtox = @(x) lunique(cums(1:floor(N*x)))/out.totnuq;
out.nuqp1 = cumtox(0.01);
out.nuqp10 = cumtox(0.1);
out.nuqp20 = cumtox(0.2);
out.nuqp50 = cumtox(0.5);

% (**1**) how many unique extrema are in first <length> of time series

Ns = [10, 50, 100, 1000];
for i = 1:length(Ns)
    if N >= Ns(i)
        out.(sprintf('nuql%u',Ns(i))) = lunique(cums(1:Ns(i)))/out.totnuq;
    else
        out.(sprintf('nuql%u',Ns(i))) = NaN;
    end
end

% (**2**) Actual proportion of full range captured at different points

out.p1 = cums(ceil(N*0.01))/fullr;
out.p10 = cums(ceil(N*0.1))/fullr;
out.p20 = cums(ceil(N*0.2))/fullr;
out.p50 = cums(ceil(N*0.5))/fullr;


Ns = [10, 50, 100, 1000];
for i = 1:length(Ns)
    if N >= Ns(i)
        out.(sprintf('l%u',Ns(i))) = cums(Ns(i))/fullr;
    else
        out.(sprintf('l%u',Ns(i))) = NaN;
    end
end

end
