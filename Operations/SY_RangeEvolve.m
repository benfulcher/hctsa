% SY_RangeEvolve
% 
% Measures of the range of the time series as a function of time,
% i.e., range(x_{1:i}) for i = 1, 2, ..., N, where N is the length of the time
% series.
% 
% INPUT: y, the time series
% 
% Outputs are based on the dynamics of how new extreme events occur with time.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = SY_RangeEvolve(y)
% Ben Fulcher, September 2009

doplot = 0; % plot outputs
N = length(y); % length of the time series
cums = zeros(N,1);

for i = 1:N
    cums(i) = range(y(1:i));
end
% cums=cums/range(y);

if doplot
    figure('color','w');
    plot(cums);
end

fullr = range(y);

% Define an anonymous function for the number of unique entries in a vector, x:
lunique = @(x) length(unique(x));

out.totnuq = lunique(cums);

% how many of the unique extrema are in first <proportion> of time series
cumtox = @(x) lunique(cums(1:floor(N*x)))/out.totnuq;
out.nuqp1 = cumtox(0.01);
out.nuqp10 = cumtox(0.1);
out.nuqp20 = cumtox(0.2);
out.nuqp50 = cumtox(0.5);
% out.nuqp1 = lunique(cums(1:floor(N*0.01)))/out.totnuq;
% out.nuqp10 = lunique(cums(1:floor(N*0.1)))/out.totnuq;
% out.nuqp20 = lunique(cums(1:floor(N*0.2)))/out.totnuq;
% out.nuqp50 = lunique(cums(1:floor(N*0.5)))/out.totnuq;

% how many unique extrema are in first <length> of time series
Ns = [10, 50, 100, 1000];
for i = 1:length(Ns)
    if N >= Ns(i)
        herehere = lunique(cums(1:Ns(i)))/out.totnuq;
        eval(sprintf('out.nuql%u = herehere;',Ns(i)));
    else
        eval(sprintf('out.nuql%u = NaN;',Ns(i)));
    end
end
% if N > 10
%     out.nuql10 = lunique(cums(1:10))/out.totnuq;
% else
%     out.nuql10 = NaN;
% end
% if N > 50
%     
% out.nuql50 = lunique(cums(1:50))/out.totnuq;
% 
% if N > 100
%     out.nuql100 = lunique(cums(1:100))/out.totnuq;
% else
%     out.nuql100 = NaN;
% end
%     
% if N > 1000
%     out.nuql1000 = lunique(cums(1:1000))/out.totnuq;
% else
%     out.nuql1000 = NaN;
% end


% (**2**) Actual proportion of full range captured at different points

out.p1 = cums(ceil(N*0.01))/fullr;
out.p10 = cums(ceil(N*0.1))/fullr;
out.p20 = cums(ceil(N*0.2))/fullr;
out.p50 = cums(ceil(N*0.5))/fullr;


Ns = [10, 50, 100, 1000];
for i = 1:length(Ns)
    if N >= Ns(i)
        herehere = cums(Ns(i))/fullr;
        eval(sprintf('out.l%u = herehere;',Ns(i)));
    else
        eval(sprintf('out.l%u = NaN;',Ns(i)));
    end
end
% out.l10 = cums(10)/fullr;
% out.l50 = cums(50)/fullr;
% 
% if N > 100
%     out.l100 = cums(100)/fullr;
% else
%     out.l100 = NaN;
% end
% 
% if N > 1000,
%     out.l1000 = cums(1000)/fullr;
% else
%     out.l1000 = NaN;
% end

end