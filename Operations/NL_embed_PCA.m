function out = NL_embed_PCA(y,tau,m)
% NL_embed_PCA  Principal Components analysis of a time series in an embedding space.
%
% Reconstructs the time series as a time-delay embedding, and performs Principal
% Components Analysis on the result using princomp code from
% Matlab's Bioinformatics Toolbox.
%
% This technique is known as singular spectrum analysis.
%
% "Extracting qualitative dynamics from experimental data"
% D. S. Broomhead and G. P. King, Physica D 20(2-3) 217 (1986)
%
%---INPUTS:
% y, the input time series
%
% tau, the time-delay, can be an integer or 'ac', or 'mi' for first
%               zero-crossing of the autocorrelation function or first minimum
%               of the automutual information, respectively
%
% m, the embedding dimension
%
% OUTPUTS: Various statistics summarizing the obtained eigenvalue distribution.
%
% The suggestion to implement this idea was provided by Siddarth Arora.
% (Siddharth Arora, <arora@maths.ox.ac.uk>)

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
% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(tau)
    tau = 'ac'; % embed by first zero-crossing of autocorrelation function
end

if nargin < 3 || isempty(m)
    m = 3; % three-dimensional embedding
end

% Embed the signal via time-delay method
y_embed = BF_embed(y,tau,m,0);

if isnan(y_embed);
    % embedding parameters are unsuitable (likely that tau is too long...)
    fprintf(1,'Embedding parameters are not suitable for this time series\n');
    out = NaN; return
end

doPlot = 0;

% ------------------------------------------------------------------------------
% Do the pca using Statistics toolbox function, 'princomp'
% ------------------------------------------------------------------------------
[~, ~, latent] = pca(y_embed);

perc = latent/sum(latent); % proportion of variance explained

if doPlot
    figure('color','w'); plot(perc,'o-k'); ylim([0,1]);
end

%-------------------------------------------------------------------------------
% Raw outputs:
%-------------------------------------------------------------------------------
for i = 1:m
    out.(sprintf('perc_%u',i)) = perc(i);
end

% ------------------------------------------------------------------------------
%% Get statistics of the eigenvalue distribution
% ------------------------------------------------------------------------------
out.std = std(perc);
out.range = max(perc) - min(perc);
out.min = min(perc); % (this is the same as perc_(m), since perc is decreasing)
out.max = max(perc); % (this is the same as perc_1)
out.top2 = sum(perc(1:2)); % variance explained in top two eigendirections

% Number of eigenvalues you need to reconstruct X%
csperc = cumsum(perc);
out.nto50 = first_fn(csperc,0.5,+1);
out.nto60 = first_fn(csperc,0.6,+1);
out.nto70 = first_fn(csperc,0.7,+1);
out.nto80 = first_fn(csperc,0.8,+1);
out.nto90 = first_fn(csperc,0.9,+1);

% When individual % variance explained goes below X for the first time:
out.fb05 = first_fn(perc,0.5,-1);
out.fb02 = first_fn(perc,0.2,-1);
out.fb01 = first_fn(perc,0.1,-1);
out.fb001 = first_fn(perc,0.01,-1);

% ------------------------------------------------------------------------------
function firsti = first_fn(p,threshold,overOrUnder)
    % Find the first time p goes under the threshold, x

    if overOrUnder==-1
        % Under threshold
        firsti = find(p < threshold,1,'first');
    else
        % Over threshold
        firsti = find(p > threshold,1,'first');
    end

    % If it never goes under -- saturate as m at the maximum
    % (could be NaN, but this is more interpretable/comparable)
    if isempty(firsti)
        firsti = length(p) + 1;
    end

end

end
