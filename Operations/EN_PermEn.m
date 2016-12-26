function out = EN_PermEn(y,m,tau)
% EN_PermEn     Permutation Entropy of a time series.
%
% "Permutation Entropy: A Natural Complexity Measure for Time Series"
% C. Bandt and B. Pompe, Phys. Rev. Lett. 88(17) 174102 (2002)
%
%---INPUTS:
% y, the input time series
% m, the embedding dimension (or order of the permutation entropy)
% tau, the time-delay for the embedding
%
%---OUTPUT:
% Outputs the permutation entropy and normalized version computed according to
% different implementations

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

%-------------------------------------------------------------------------------
% Check inputs and set defaults:
%-------------------------------------------------------------------------------

if nargin < 2 || isempty(m)
    m = 2; % order 2
end

if nargin < 3 || isempty(tau)
    tau = 1;
end

% ------------------------------------------------------------------------------
% Embed the signal:
% ------------------------------------------------------------------------------
x = BF_embed(y,tau,m,0);
Nx = length(x); % number of embedding vectors produced

% Generate permutations up to the embedding dimension, m:
permList = perms(1:m);
numPerms = length(permList);

% Initialize
countPerms = zeros(numPerms,1);

% Count each type of permutation through the time series
for j = 1:Nx

    % Get the permutation for this local time-series segment:
    [~,ix] = sort(x(j,:));

    % Match this to one of the permutations:

    % (i) Nicer but slower:
    % thisPerm = find(all(bsxfun(@minus,ix,permList)==0,2),1);
    % countPerms(thisPerm) = countPerms(thisPerm) + 1;

    % (ii) Uglier but faster:
    for k = 1:numPerms
        if all(permList(k,:)-ix == 0)
            % We found it! Increment:
            countPerms(k) = countPerms(k) + 1;
            break
        end
    end
end

% ------------------------------------------------------------------------------
% Convert counts to probabilities
p = countPerms/Nx; %((Nx-(m-1))*tau);
p_0 = p(p>0); % makes log(0) = 0
out.permEn = -sum(p_0.*log2(p_0));

% Normalized permutation entropy (more comparable across m?)
mFact = factorial(m);
out.normPermEn = out.permEn/log2(mFact);

% ------------------------------------------------------------------------------
% Adapted implementation by Bruce Land and Damian Elias:
% cf.:
% http://people.ece.cornell.edu/land/PROJECTS/Complexity/
% http://people.ece.cornell.edu/land/PROJECTS/Complexity/logisticPE.m

% Not clear to me why you would make log(0) = log(1/Nx); (the minimum)
% rather than exclude it from the sum, as is done here:
p_LE = max(1/Nx,p);
out.permEnLE = -sum(p_LE.*log(p_LE))/(m-1);

end
