function out = MD_pNN(x)
% MD_pNN    pNNx measures of heart rate variability.
%
% Applies pNNx measures to time series assumed to represent sequences of
% consecutive RR intervals measured in milliseconds.
%
% cf. "The pNNx files: re-examining a widely used heart rate variability
%           measure", J.E. Mietus et al., Heart 88(4) 378 (2002)
%
%---INPUTS:
% x, the input time series
%
% This code is derived from MD_hrv_classic.m becuase it doesn't make medical
% sense to do PNN on a z-scored time series.
%
% But now PSD doesn't make too much sense, so we just evaluate the pNN measures.
%
% Code is heavily derived from that provided by Max A. Little:
% http://www.maxlittle.net/

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Max A. Little, <max.a.little@gmail.com>,
% <http://www.maxlittle.net/> and Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

diffx = diff(x); % successive increments of the time series
N = length(x); % length of the time series

% ------------------------------------------------------------------------------
%% Calculate pNNx percentage
% ------------------------------------------------------------------------------
% strange to do this for a z-scored time series...
% pnntime = 20;

Dx = abs(diffx) * 1000; % assume milliseconds as for RR intervals
pnns = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
for i = 1:length(pnns)
    out.(sprintf('pnn%u',pnns(i))) = sum(Dx > pnns(i)) / (N-1);
end

end
