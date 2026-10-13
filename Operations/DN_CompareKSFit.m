function out = DN_CompareKSFit(x, whatDistn)
% DN_CompareKSFit   Compares a fitted distribution with the smoothed distribution of the data.
%
% Fits a standard distribution to the data (by maximum likelihood, using
% Statistics Toolbox functions), and compares it with a kernel-smoothed
% estimate of the distribution of the values (BF_KSDensity, a Gaussian kernel with
% a normal-reference bandwidth; "KS" here means kernel-smoothed, not
% Kolmogorov-Smirnov). Both curves are evaluated on a
% common grid of 1000 points that covers the smoothed distribution and the body
% of the fitted distribution (out to where it falls to 1/100 of its peak, found
% by root finding; or, for a positive-only distribution whose density is
% unbounded or maximal at 0, from 0 to its 0.99 quantile). The smoothed
% distribution is evaluated at the grid points, and the fitted distribution as
% its mean density over each grid cell (its probability in the cell divided by
% the cell width), which stays finite where the density diverges. They are then
% compared by the area between them, the separation of their peaks, their
% overlap, and the relative entropy.
%
% The exponential and Rayleigh distributions require non-negative values, and
% the gamma, log-normal and Weibull distributions require positive values;
% NaN is returned if the data do not satisfy this (and for a constant series
% in the Rayleigh and exponential cases). For the beta distribution, the data
% are first rescaled to lie inside (0, 1), and all outputs are then in rescaled
% units.
%
% ---INPUTS:
% x, the input data vector
% whatDistn, the type of distribution to fit to the data:
%           'norm' (Gaussian; normfit), 'ev' (extreme value; evfit),
%           'uni' (uniform; unifit), 'beta' (beta; betafit),
%           'rayleigh' (Rayleigh; raylfit), 'exp' (exponential; expfit),
%           'gamma' (gamma; gamfit), 'logn' (log-normal; lognfit),
%           'wbl' (Weibull; wblfit).
%
% ---OUTPUTS:
% adiff, the absolute area between the two distributions (0 for a perfect
%       match, at most 2).
% peaksepy, the separation in height between the peaks: the maximum of the
%       fitted distribution minus that of the smoothed distribution.
% peaksepx, the separation in position between the peaks: the position of the
%       peak of the fitted distribution minus that of the smoothed distribution
%       (NaN for the uniform, which has no peak).
% olapint, the overlap integral of the two distributions, multiplied by the
%       standard deviation of the data so that it does not depend on their
%       scale.
% relent, the relative entropy (Kullback-Leibler divergence), in nats, of the
%       fitted distribution from the smoothed distribution.
%
% ---NOTES:
% adiff, olapint and relent do not depend on the scale of the data, but
% peaksepy and peaksepx do.

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

% ------------------------------------------------------------------------------
%% Preprocessing
% ------------------------------------------------------------------------------

% Fit the distribution 'whatDistn' to the input data, x.
xStep = std(x) / 100; % initial stride when searching outward for the edges of the fitted distribution
q99 = []; % right edge set by a quantile (densities unbounded or maximal at zero)

switch whatDistn
	case 'norm'
		[a, b] = normfit(x);
		pdfFn = @(z) normpdf(z, a, b);
		cdfLo = @(z) normcdf(z, a, b);
		cdfUp = @(z) normcdf(z, a, b, 'upper');
		thresh = pdfFn(a) / 100; % the edges are where the density falls to 1/100 of its peak
		leftStart = a; rightStart = a;

	case 'ev'
		a = evfit(x);
		pdfFn = @(z) evpdf(z, a(1), a(2));
		cdfLo = @(z) evcdf(z, a(1), a(2));
		cdfUp = @(z) evcdf(z, a(1), a(2), 'upper');
		thresh = pdfFn(a(1)) / 100;
		leftStart = a(1); rightStart = a(1);

	case 'uni'
		[a, b] = unifit(x);
		pdfFn = @(z) unifpdf(z, a, b);
		cdfLo = @(z) unifcdf(z, a, b);
		cdfUp = @(z) unifcdf(z, a, b, 'upper');
		thresh = pdfFn(mean(x)) / 100;
		leftStart = mean(x); rightStart = mean(x);

	case 'beta'
		% clumsily scale to the range (0,1)
		x = (x - min(x) + 0.01 * std(x)) / (max(x) - min(x) + 0.02 * std(x));
		xStep = std(x) / 100; % will need a new step size for the rescaled data
		a = betafit(x);
		pdfFn = @(z) betapdf(z, a(1), a(2));
		cdfLo = @(z) betacdf(z, a(1), a(2));
		cdfUp = @(z) betacdf(z, a(1), a(2), 'upper');
		thresh = 1E-5; % ok -- consistent since all scaled to the same range
		leftStart = mean(x); rightStart = mean(x);

	case 'rayleigh'
		if any(x < 0)
			fprintf(1, 'The data are not positive, but Rayleigh is a positive-only distribution.\n');
			out = NaN; return
		elseif all(x == x(1))
			fprintf(1, 'Data are a constant\n');
			out = NaN; return
		end
		a = raylfit(x);
		pdfFn = @(z) raylpdf(z, a);
		cdfLo = @(z) raylcdf(z, a);
		cdfUp = @(z) raylcdf(z, a, 'upper');
		thresh = pdfFn(a) / 100; % (the peak is at the scale parameter)
		leftStart = []; rightStart = a;

	case 'exp'
		if any(x < 0)
			fprintf(1, 'The data contains negative values, but Exponential is a positive-only distribution.\n');
			out = NaN; return
		elseif all(x == x(1))
			fprintf(1, 'Data are a constant\n');
			out = NaN; return
		end
		a = expfit(x);
		pdfFn = @(z) exppdf(z, a);
		cdfLo = @(z) expcdf(z, a);
		cdfUp = @(z) expcdf(z, a, 'upper');
		thresh = pdfFn(0) / 100; % (the peak is at zero)
		leftStart = []; rightStart = 0;

	case 'gamma'
		if any(x <= 0)
			% (the gamma likelihood is degenerate at an exact zero: log(0) enters it)
			fprintf(1, 'The data are not positive, but Gamma is a positive-only distribution.\n');
			out = NaN; return
		end
		a = gamfit(x);
		pdfFn = @(z) gampdf(z, a(1), a(2));
		cdfLo = @(z) gamcdf(z, a(1), a(2));
		cdfUp = @(z) gamcdf(z, a(1), a(2), 'upper');
		leftStart = [];
		if a(1) < 1 % the density is unbounded at zero
			q99 = gaminv(0.99, a(1), a(2));
		else
			rightStart = (a(1) - 1) * a(2); % the mode
			thresh = pdfFn(rightStart) / 100;
		end

	case 'logn'
		if any(x <= 0)
			fprintf(1, 'The data are not positive, but Log-Normal is a positive-only distribution.\n');
			out = NaN; return
		end
		a = lognfit(x);
		pdfFn = @(z) lognpdf(z, a(1), a(2));
		cdfLo = @(z) logncdf(z, a(1), a(2));
		cdfUp = @(z) logncdf(z, a(1), a(2), 'upper');
		rightStart = exp(a(1) - a(2)^2); % the mode
		thresh = pdfFn(rightStart) / 100;
		leftStart = [];

	case 'wbl'
		if any(x <= 0)
			fprintf(1, 'The data are not positive, but Weibull is a positive-only distribution.\n');
			out = NaN; return
		end
		a = wblfit(x);
		pdfFn = @(z) wblpdf(z, a(1), a(2));
		cdfLo = @(z) wblcdf(z, a(1), a(2));
		cdfUp = @(z) wblcdf(z, a(1), a(2), 'upper');
		leftStart = [];
		if a(2) <= 1 % the density is unbounded (shape < 1) or maximal (shape = 1) at zero
			q99 = wblinv(0.99, a(1), a(2));
		else
			rightStart = a(1) * ((a(2) - 1) / a(2))^(1 / a(2)); % the mode
			thresh = pdfFn(rightStart) / 100;
		end

	otherwise
		error('Unknown distribution: %s.', whatDistn)
end

% The range [xf(1), xf(2)] of the body of the fitted distribution: where its density
% falls to thresh, found by root finding outward from the peak; for the positive-only
% distributions the left edge is 0, and the right edge is the 0.99 quantile when the
% density is unbounded or maximal at 0 (for an exponential, exactly where it falls to
% 1/100 of its peak).
if ~isempty(q99)
	xf = [0, q99];
elseif isempty(leftStart)
	xf = [0, SUB_DensityEdge(pdfFn, rightStart, 1, xStep, thresh)];
else
	xf = [SUB_DensityEdge(pdfFn, leftStart, -1, xStep, thresh), ...
			SUB_DensityEdge(pdfFn, rightStart, 1, xStep, thresh)];
end

% ------------------------------------------------------------------------------
% Estimate smoothed empirical distribution
% ------------------------------------------------------------------------------
[f, xi] = BF_KSDensity(x);
xi = xi(f > 1E-6); % only keep values greater than 1E-6
if isempty(xi)
	out = NaN; return
	% in future could change to the threshold from 1E-6 to some fraction
	% of the peak value...
end

xi = [floor(xi(1) * 10) / 10, ceil(xi(end) * 10) / 10];

% Find appropriate range [x1 x2] that incorporates the full range of both
x1 = min([xf(1), xi(1)]);
x2 = max([xf(2), xi(end)]);

% ------------------------------------------------------------------------------
%% Evaluate both over the same grid
% ------------------------------------------------------------------------------
xi = linspace(x1, x2, 1000);
dx = xi(2) - xi(1);
f = BF_KSDensity(x, xi); % the smoothed empirical distribution
% The fitted distribution as its mean density over each grid cell [xi - dx/2, xi + dx/2]
% (its probability in the cell divided by dx), which stays finite where the density
% diverges (a gamma or Weibull with shape < 1 at 0) and does not depend on where the
% grid points fall relative to a jump (e.g., the ends of a uniform). The probability
% is a difference of the lower-tail CDF below the median and of the upper-tail CDF
% above it, so that tail cells keep their relative precision.
lower = cdfLo(xi + dx/2) - cdfLo(xi - dx/2);
upper = cdfUp(xi - dx/2) - cdfUp(xi + dx/2);
useLower = cdfLo(xi) <= 0.5;
ffit = zeros(size(xi));
ffit(useLower) = lower(useLower) / dx;
ffit(~useLower) = upper(~useLower) / dx;
ffit = max(ffit, 0); % (rounding)

% now the two cover the same range in x

% ------------------------------------------------------------------------------
%% Outputs:
% ------------------------------------------------------------------------------

% ADIFF: returns absolute area between the curves
out.adiff = sum(abs(f - ffit) * dx);

% PEAKSEPY: returns the seperation (in y) between the maxima of each distrn
% NOT A GREAT STATISTIC FOR GENERAL DISTRIBUTIONS -- scales with variance
max1 = max(f);
max2 = max(ffit);
out.peaksepy = max2 - max1;

% PEAKSEPX: returns the seperation (in x) between the maxima of each distrn
% NOT A VERY WELL POSED STATISTIC FOR GENERAL DISTRIBUTIONS since this scales with variance
if strcmp(whatDistn, 'uni')
	out.peaksepx = NaN; % a uniform density is flat: it has no peak position
else
	[~, i1] = max(f);
	[~, i2] = max(ffit);
	out.peaksepx = xi(i2) - xi(i1);
end

% OLAPINT: returns the overlap integral between the two curves; multiplying by
% std(x) makes this scale-invariant (a PDF's overlap integral scales as 1/std
% under rescaling of x, so this cancels that out -- multiplying by var(x)
% instead would not be scale-invariant)
out.olapint = sum(f .* ffit * dx) * std(x);

% RELENT: returns the relative entropy of the two distributions
% By convention 0*log(0/y) := 0, so also exclude points where f == 0 -- otherwise
% log(0/ffit) = -Inf and 0*(-Inf) evaluates to NaN, which propagated into a real fraction
% of series (up to ~9% for some distributions) rather than contributing zero as intended.
r = (ffit ~= 0) & (f ~= 0);
out.relent = sum(f(r) .* log(f(r) ./ ffit(r)) * dx);

end

% ------------------------------------------------------------------------------
function xEdge = SUB_DensityEdge(pdfFn, xStart, direction, xStep, thresh)
	% The point beyond xStart (direction = +1 or -1) at which a unimodal density falls to
	% thresh: the distance from xStart is doubled, from xStep, until the density is at most
	% thresh, and the crossing in that last interval is found by root finding (fzero; for a
	% density with a jump, such as the uniform, this is the jump). xStart if the density
	% is already at most thresh there.
	if ~(pdfFn(xStart) > thresh)
		xEdge = xStart;
		return
	end
	inside = xStart;
	stride = xStep;
	for k = 1:200
		outside = xStart + direction * stride;
		if ~isfinite(outside)
			break
		end
		if ~(pdfFn(outside) > thresh)
			xEdge = fzero(@(z) pdfFn(z) - thresh, sort([inside, outside]));
			return
		end
		inside = outside;
		stride = 2 * stride;
	end
	error('Could not bracket the edge of the fitted distribution');
end
