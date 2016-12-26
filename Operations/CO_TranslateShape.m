function out = CO_TranslateShape(y,shape,d,howToMove)
% CO_TranslateShape     Statistics on the number of datapoints residing inside
% geometric shapes moved across the time series.
%
% Inputs specify a shape and its size, and a method for moving this shape
% through the time domain.
%
% This is usually more informative in an embedding space (CO_Embed2_...), but
% here we do it just in the temporal domain (_t_).
%
% In the future, could perform a similar analysis with a soft boundary, some
% decaying force function V(r), or perhaps truncated...?
%
%---INPUTS:
% y, the input time series
% shape, the shape to move about the time-domain ('circle')
% d, a parameter specifying the size of the shape (e.g., d = 2)
% howToMove, a method specifying how to move the shape about, e.g., 'pts'
%               places the shape on each point in the time series.

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
%% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(shape)
    shape = 'circle';
end
if nargin < 3 || isempty(d)
    d = 2; % a default distance d = 2
end
if nargin < 4 || isempty(howToMove)
    howToMove = 'pts'; % by default, places shapes on each timepoint
end

% ------------------------------------------------------------------------------
% Preliminaries:
% ------------------------------------------------------------------------------
N = length(y); % the length of the time series

% y must be a column vector, transpose it if it's a row vector
if size(y,2) > size(y,1)
    y = y';
end

% Add a time index
ty = [(1:N)', y]; % has increasing integers as time in the first column

%-------------------------------------------------------------------------------
% Generate the statistics on the number of points inside the shape as it is
% translated across the time series
%-------------------------------------------------------------------------------
switch howToMove
    case 'pts' % Place shapes on each timepoint (excluding a range at start and end)
        switch shape
            case 'circle' % uses a circle of radius 'd'
                r = d;
                w = floor(r); % only consider a window radius w (these are the
                            %    only points that could possibly be inside)
                rnge = 1+w:N-w;
                NN = length(rnge); % number of admissible points
                np = zeros(NN,1); % number of points
                for i = 1:NN
                    win = ty(rnge(i)-w:rnge(i)+w,:); % create window
                    difwin = win - ones(2*w+1,1)*ty(rnge(i),:);
                    np(i) = sum(sum(difwin.^2,2) <= r^2); % number of points enclosed in shape
                end
            case 'rectangle'
                % uses a rectangle of half-width d (integer), height as double the current time
                % series value (on either side of origin)
                w = d;
                rnge = 1+d:N-d;
                NN = length(rnge); % number of admissible points
                np = zeros(NN,1); % number of points
                for i = 1:NN
                    np(i) = sum(abs(y(rnge(i)-d:rnge(i)+d)) <= abs(y(i)));
                end
        otherwise
            error('Unknown shape ''%s''',shape)
        end
    otherwise
        error('Unknown setting for ''howToMove'' input: ''%s''',howToMove)
end

%-------------------------------------------------------------------------------

% -----
% Maximum possible hits in the shape:
% -----
out.max = max(np);
out.std = std(np);
out.mean = mean(np);

% Count the hits:
histnp = arrayfun(@(x)sum(np==x),unique(np));
% bar(histnp)

% Compute mode:
[out.npatmode, out.mode] = max(histnp);
out.npatmode = out.npatmode/NN;

% Output all stats:
if 2*w + 1 >= 1; out.ones = mean(np==1); end
if 2*w + 1 >= 2; out.twos = mean(np==2); end
if 2*w + 1 >= 3; out.threes = mean(np==3); end
if 2*w + 1 >= 4; out.fours = mean(np==4); end
if 2*w + 1 >= 5; out.fives = mean(np==5); end
if 2*w + 1 >= 6; out.sixes = mean(np==6); end
if 2*w + 1 >= 7; out.sevens = mean(np==7); end
if 2*w + 1 >= 8; out.eights = mean(np==8); end
if 2*w + 1 >= 9; out.nines = mean(np==9); end
if 2*w + 1 >= 10; out.tens = mean(np==10); end
if 2*w + 1 >= 11; out.elevens = mean(np==11); end

% -----
% Stationarity in 2,3,4 segments
% -----
out.statav2_m = SY_SlidingWindow(np,'mean','std',2,1);
out.statav2_s = SY_SlidingWindow(np,'std','std',2,1);
out.statav3_m = SY_SlidingWindow(np,'mean','std',3,1);
out.statav3_s = SY_SlidingWindow(np,'std','std',3,1);
out.statav4_m = SY_SlidingWindow(np,'mean','std',4,1);
out.statav4_s = SY_SlidingWindow(np,'std','std',4,1);

% plot(np)



end
