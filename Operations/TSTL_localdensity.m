function out = TSTL_localdensity(y,NNR,past,embedParams)
% TSTL_localdensity     Local density estimates in the time-delay embedding space
%
% TSTOOL code localdensity is very poorly documented in the TSTOOL
% package, can assume it returns local density estimates in the
% time-delay embedding space.
%
% TSTOOL: http://www.physik3.gwdg.de/tstool/
%
%---INPUTS:
%
% y, the time series as a column vector
%
% NNR, number of nearest neighbours to compute
%
% past, number of time-correlated points to discard (samples)
%
% embedParams, the embedding parameters, inputs to BF_embed as {tau,m}, where
%               tau and m can be characters specifying a given automatic method
%               of determining tau and/or m (see BF_embed).
%
%---OUTPUTS: various statistics on the local density estimates at each point in
% the time-delay embedding, including the minimum and maximum values, the range,
% the standard deviation, mean, median, and autocorrelation.

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
%% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(NNR)
    NNR = 3; % 3 nearest neighbours
end

if nargin < 3 || isempty(past)
    past = 40;
end

if nargin < 4 || isempty(embedParams)
    embedParams = {'ac','fnnmar'};
    fprintf(1,'Using default embedding using autocorrelation and cao''s method.\n');
end

% ------------------------------------------------------------------------------
%% Embed the signal
% ------------------------------------------------------------------------------
s = BF_embed(y,embedParams{1},embedParams{2},1);

if ~isa(s,'signal') && isnan(s); % embedding failed
    error('Embedding failed.')
end

% ------------------------------------------------------------------------------
%% Run the TSTOOL function, localdensity
% ------------------------------------------------------------------------------
rs = localdensity(s,NNR,past);

% ------------------------------------------------------------------------------
%% Convert output to data
% ------------------------------------------------------------------------------
locden = data(rs);
if all(locden == 0)
    out = NaN; return
end
% locden is a vector of length equal to the number of points in the
% embedding space (length of time series - m + 1), presumably the local
% at each point

out.minden = min(locden);
out.maxden = max(locden);
out.iqrden = iqr(locden);
out.rangeden = range(locden);
out.stdden = std(locden);
out.meanden = mean(locden);
out.medianden = median(locden);

F_acden = @(x) CO_AutoCorr(locden,x,'Fourier'); % autocorrelation of locden for 1:5
for i = 1:5
    out.(sprintf('ac%uden',i)) = F_acden(i);
end

% Estimates of correlation length:
out.tauacden = CO_FirstZero(locden,'ac'); % first zero-crossing of autocorrelation function
out.taumiden = CO_FirstMin(locden,'mi'); % first minimum of automutual information function

end
