% NL_embed_PCA
% 
% Reconstructs the time series as a time-delay embedding, and performs Principal
% Components Analysis on the result using princomp code from
% Matlab's Bioinformatics Toolbox.
% 
% This technique is known as singular spectrum analysis
% 
% "Extracting qualitative dynamics from experimental data"
% D. S. Broomhead and G. P. King, Physica D 20(2-3) 217 (1986)
% 
% INPUTS:
% y, the input time series
% 
% tau, the time-delay, can be an integer or 'ac', or 'mi' for first
%               zero-crossing of the autocorrelation function or first minimum
%               of the automutual information, respectively
%               
% m, the embedding dimension
% 
% Outputs are various statistics summarizing the obtained eigenvalue distribution.
% 
% The suggestion to implement this idea was provided by Siddarth Arora.
% (Siddharth Arora, <arora@maths.ox.ac.uk>)
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

function out = NL_embed_PCA(y,tau,m)
% Ben Fulcher, 25/2/2010

if nargin < 2 || isempty(tau)
    tau = 'ac'; % embed by first zero-crossing of autocorrelation function
end

if nargin < 3 || isempty(m)
    m = 3; % three-dimensional embedding
end

% embed the signal via time-delay method
y_embed = BF_embed(y,tau,m,0);

if isnan(y_embed);
    % embedding parameters are unsuitable (likely that tau is too long...)
    fprintf(1,'Embedding parameters are not suitable for this time series\n');
    out = NaN; return
end

% Do the pca using Bioinformatics toolbox routine 'princomp'
[pc, score, latent] = princomp(y_embed);

% perc=round(latent/sum(latent)*1000)/10; % percentage of variance explained (1 d.p.)
% sum(latent)
perc = latent/sum(latent); % percentage of variance explained

% plot(perc); ylim([0,1]);

%% Get statistics off the eigenvalue distribution
csperc = cumsum(perc);
out.std = std(perc);
out.range = max(perc) - min(perc);
out.max = max(perc);
out.min = min(perc);
out.top2 = sum(perc(1:2)); % variance explained in top two eigendirections
out.nto80 = find(csperc>0.8,1,'first'); % number of eigenvalues you need to reconstruct 80%
out.nto90 = find(csperc>0.9,1,'first'); % number of eigenvalues you need to reconstruct 90%

out.fb01 = find(perc < 0.1,1,'first'); % when perc goes below 0.1 for the first time
if isempty(out.fb01), out.fb01 = length(perc) + 1; end % could make it NaN...

out.fb001 = find(perc < 0.01,1,'first'); % when perc goes below 0.01 for the first time
if isempty(out.fb001), out.fb001 = length(perc)+1; end % could also make it NaN...


end