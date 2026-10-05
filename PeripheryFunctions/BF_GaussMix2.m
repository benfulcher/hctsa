function [w, mu, sigma] = BF_GaussMix2(xc, p, delta)
% BF_GaussMix2   Mixture of two Gaussians fitted to a binned distribution, by EM.
%
% Fits the density w(1)*N(mu(1), sigma(1)^2) + w(2)*N(mu(2), sigma(2)^2) to a
% distribution given as a histogram or density on an equally spaced grid, by the
% expectation-maximization algorithm applied to the grid points, weighted by the
% mass p at each point (maximum likelihood for the binned data). The variance of
% each component includes delta^2/12, the variance within a bin of width delta, so
% that a component cannot be narrower than a bin (as for a point mass at a bin
% center) and the fitted density is consistent with the resolution of the grid.
% The start is deterministic (no random initialization): the two means start at
% the lower and upper quartiles of the distribution, the weights at 1/2 and the
% standard deviations at that of the distribution. The algorithm runs until the
% log-likelihood changes by less than 1e-10 (at most 1000 iterations). Where the
% distribution is not bimodal the two components are not identifiable (many
% mixtures describe the density equally well), but the fitted density is.
%
% ---INPUTS:
% xc, the grid points (e.g., bin centers), a column vector
% p, the density (or counts) at each point, a column vector of non-negative values
%       (rescaled here to sum to 1)
% delta, the grid spacing (bin width)
%
% ---OUTPUTS: 1 x 2 vectors, with the components ordered by increasing mean
% w, the mixture weights
% mu, the means
% sigma, the standard deviations

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
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

xc = xc(:);
p = p(:) / sum(p); % mass at each grid point
v0 = delta^2/12; % variance within a bin

% Start: means at the quartiles, equal weights, the overall standard deviation
cp = cumsum(p);
mu = [xc(find(cp >= 0.25, 1)), xc(find(cp >= 0.75, 1))];
mBar = sum(p.*xc);
sBar = sqrt(sum(p.*(xc - mBar).^2) + v0);
if mu(1) == mu(2) % both quartiles in one bin
	mu = mBar + sBar*[-0.5, 0.5];
end
w = [0.5, 0.5];
sigma = sBar*[1, 1];

llOld = -Inf;
for iter = 1:1000
	% E step: responsibility of component 1 at each grid point (and the log-likelihood)
	l1 = log(w(1)) - log(sigma(1)) - (xc - mu(1)).^2 / (2*sigma(1)^2);
	l2 = log(w(2)) - log(sigma(2)) - (xc - mu(2)).^2 / (2*sigma(2)^2);
	lmax = max(l1, l2);
	ll = sum(p.*(lmax + log(exp(l1 - lmax) + exp(l2 - lmax))));
	r1 = 1 ./ (1 + exp(l2 - l1));
	if abs(ll - llOld) < 1e-10, break; end
	llOld = ll;

	% M step
	n1 = sum(p.*r1); n2 = 1 - n1;
	if n1 < 1e-8 || n2 < 1e-8, break; end % one component has vanished
	w = [n1, n2];
	mu = [sum(p.*r1.*xc)/n1, sum(p.*(1 - r1).*xc)/n2];
	sigma = sqrt([sum(p.*r1.*(xc - mu(1)).^2)/n1, sum(p.*(1 - r1).*(xc - mu(2)).^2)/n2] + v0);
end

[mu, ix] = sort(mu); % order the components by mean
w = w(ix);
sigma = sigma(ix);

end
