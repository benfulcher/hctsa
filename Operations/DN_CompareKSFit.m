function out = DN_CompareKSFit(x,whatDistn)
% DN_CompareKSFit       Fits a distribution to data.
%
% Returns simple statistics on the discrepancy between the
% kernel-smoothed distribution of the time-series values, and the distribution
% fitted to it by some model: Gaussian (using normfifit from Matlab's
% Statistics Toolbox), Extreme Value (evfifit), Uniform (unififit), Beta
% (betafifit), Rayleigh (raylfifit), Exponential (expfifit), Gamma (gamfit),
% LogNormal (lognfifit), and Weibull (wblfifit).
%
%---INPUTS:
% x, the input data vector
% whatDistn, the type of distribution to fit to the data:
%           'norm' (normal), 'ev' (extreme value), 'uni' (uniform),
%           'beta' (Beta), 'rayleigh' (Rayleigh), 'exp' (exponential),
%           'gamma' (Gamma), 'logn' (Log-Normal), 'wbl' (Weibull).
%
%---OUTPUTS: include the absolute area between the two distributions, the peak
% separation, overlap integral, and relative entropy.

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
%% Preprocessing
% ------------------------------------------------------------------------------

% Fit the distribution 'whatDistn' to the input data, x.
xStep = std(x)/100; % set a step size

switch whatDistn
    case 'norm'
        [a, b] = normfit(x);
		peaky = normpdf(a,a,b); thresh = peaky/100; % stop when gets to 1/100 of peak value
		xf(1) = mean(x); ange = 10;
        while ange > thresh, xf(1) = xf(1)-xStep; ange = normpdf(xf(1),a,b); end
		xf(2) = mean(x); ange = 10;
        while ange > thresh, xf(2) = xf(2)+xStep; ange = normpdf(xf(2),a,b); end

    case 'ev'
        a = evfit(x);
		peaky = evpdf(a(1),a(1),a(2)); thresh = peaky/100;
		xf(1) = 0; ange = 10;
        while ange > thresh, xf(1) = xf(1)-xStep; ange = evpdf(xf(1),a(1),a(2)); end
		xf(2) = 0; ange = 10;
        while ange > thresh, xf(2) = xf(2)+xStep; ange = evpdf(xf(2),a(1),a(2)); end

    case 'uni'
        [a, b] = unifit(x);
		peaky = unifpdf(mean(x),a,b); thresh = peaky/100;
		xf(1) = 0; ange = 10;
        while ange > thresh, xf(1) = xf(1)-xStep; ange = unifpdf(xf(1),a,b); end
		xf(2) = 0; ange = 10;
        while ange > thresh, xf(2) = xf(2)+xStep; ange = unifpdf(xf(2),a,b); end

    case 'beta'
        % clumsily scale to the range (0,1)
        x = (x-min(x)+0.01*std(x))/(max(x)-min(x)+0.02*std(x));
        xStep = std(x)/100; % will need a new step size for the rescaled data
        a = betafit(x);
		thresh = 1E-5; % ok -- consistent since all scaled to the same range
		xf(1) = mean(x); ange = 10;
        while ange > thresh, xf(1) = xf(1)-xStep; ange = betapdf(xf(1),a(1),a(2)); end
		xf(2) = mean(x); ange = 10;
        while ange > thresh, xf(2) = xf(2)+xStep; ange = betapdf(xf(2),a(1),a(2)); end

    case 'rayleigh'
        if any(x < 0),
            fprintf(1,'The data are not positive, but Rayleigh is a positive-only distribution.\n');
            out = NaN;
            return
        else % fit a Rayleigh distribution to the positive-only data
            a = raylfit(x);
			peaky = raylpdf(a,a); thresh=peaky/100;
			xf(1) = 0;
			xf(2) = a; ange = 10;
            while ange > thresh, xf(2) = xf(2)+xStep; ange = raylpdf(xf(2),a); end
        end

    case 'exp'
        if any(x < 0)
            fprintf(1,'The data contains negative values, but Exponential is a positive-only distribution.\n');
            out = NaN; return
        else a = expfit(x);
			peaky = exppdf(0,a); thresh = peaky/100;
			xf(1) = 0;
			xf(2) = 0; ange = 10;
            while ange > thresh, xf(2) = xf(2)+xStep; ange = exppdf(xf(2),a); end
        end

    case 'gamma'
        if any(x < 0)
            fprintf(1,'The data contains negative values, but Gamma is a positive-only distribution.\n');
            out = NaN; return
        else a = gamfit(x);
			if a(1) < 1
				thresh = gampdf(0,a(1),a(2))/100;
			else
				peaky = gampdf((a(1)-1)*a(2),a(1),a(2)); thresh = peaky/100;
			end
			xf(1) = 0;
			xf(2) = a(1)*a(2); ange = 10;
            while ange > thresh, xf(2) = xf(2)+xStep; ange = gampdf(xf(2),a(1),a(2)); end
        end

    case 'logn'
        if any(x <= 0)
            fprintf(1,'The data are not positive, but Log-Normal is a positive-only distribution.\n');
            out = NaN; return
        else
			a = lognfit(x);
			peaky = lognpdf(exp(a(1)-a(2)^2),a(1),a(2)); thresh = peaky/100;
			xf(1) = 0;
			xf(2) = exp(a(1)-a(2)^2); ange = 10;
            while ange > thresh, xf(2) = xf(2)+xStep; ange = lognpdf(xf(2),a(1),a(2)); end
        end

    case 'wbl'
        if any(x <= 0)
            fprintf(1,'The data are not positive, but Weibull is a positive-only distribution.\n');
            out = NaN; return
        else
			a = wblfit(x);
			if a(2) <= 1;
				thresh = wblpdf(0,a(1),a(2));
			else
				peaky = wblpdf(a(1)*((a(2)-1)/a(2))^(1/a(2)),a(1),a(2));
				thresh = peaky/100;
			end
			xf(1) = 0;
			xf(2) = 0; ange = 10;
            while ange > thresh, xf(2) = xf(2)+xStep; ange = wblpdf(xf(2),a(1),a(2)); end
        end

    otherwise
        error('Unknown distribution: %s.',whatDistn)
end
% xtmafit=[floor(xtmafit(1)*10)/10 ceil(xtmafit(end)*10)/10];

% ------------------------------------------------------------------------------
% Estimate smoothed empirical distribution
% ------------------------------------------------------------------------------
[f, xi] = ksdensity(x);
xi = xi(f > 1E-6); % only keep values greater than 1E-6
if isempty(xi)
    out = NaN; return
    % in future could change to the threshold from 1E-6 to some fraction
    % of the peak value...
end

xi = [floor(xi(1)*10)/10, ceil(xi(end)*10)/10];

% Find appropriate range [x1 x2] that incorporates the full range of both
x1 = min([xf(1), xi(1)]);
x2 = max([xf(2), xi(end)]);

% ------------------------------------------------------------------------------
%% Rerun both over the same range
% ------------------------------------------------------------------------------
% (Inefficient, but easier to do it this way)
xi = linspace(x1,x2,1000);
f = ksdensity(x,xi); % the smoothed empirical distribution
switch whatDistn
    case 'norm'
        ffit = normpdf(xi,a,b);
    case 'ev'
        ffit = evpdf(xi,a(1),a(2));
    case 'uni'
        ffit = unifpdf(xi,a,b);
    case 'beta'
%         if x1<0,x1=1E-5;end
%         if x2>1,x2=1-1E-5;end
        ffit = betapdf(xi,a(1),a(2));
    case 'rayleigh'
%         if x1<0,x1=0;end
        ffit = raylpdf(xi,a);
    case 'exp'
%         if x1<0,x1=0;end
        ffit = exppdf(xi,a);
    case 'gamma'
%         if x1<0,x1=0;end
        ffit = gampdf(xi,a(1),a(2));
    case 'logn'
%         if x1<0,x1=0;end
        ffit = lognpdf(xi,a(1),a(2));
    case 'wbl'
%         if x1<0,x1=0;end
        ffit = wblpdf(xi,a(1),a(2));
end

% now the two cover the same range in x

% ------------------------------------------------------------------------------
%% Outputs:
% ------------------------------------------------------------------------------

% ADIFF: returns absolute area between the curves
out.adiff = sum(abs(f-ffit)*(xi(2)-xi(1)));

% PEAKSEPY: returns the seperation (in y) between the maxima of each distrn
% NOT A GREAT STATISTIC FOR GENERAL DISTRIBUTIONS -- scales with variance
max1 = max(f);
max2 = max(ffit);
out.peaksepy = max2-max1;

% PEAKSEPX: returns the seperation (in x) between the maxima of each distrn
% NOT A VERY WELL POSED STATISTIC FOR GENERAL DISTRIBUTIONS since this scales with variance
[~, i1] = max(f);
[~, i2] = max(ffit);
out.peaksepx = xi(i2)-xi(i1);

% OLAPINT: returns the overlap integral between the two curves; normalized by variance
out.olapint = sum(f.*ffit*(xi(2)-xi(1)))*std(x);

% RELENT: returns the relative entropy of the two distributions
r = (ffit ~= 0);
out.relent = sum(f(r).*log(f(r)./ffit(r))*(xi(2)-xi(1)));


end
