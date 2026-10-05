function out = DN_FitKernelSmooth(x, varargin)
% DN_FitKernelSmooth   Statistics of a kernel-smoothed distribution of the data.
%
% Estimates the distribution of the values with a kernel-smoothed density
% (BF_KSDensity: a Gaussian kernel with a normal-reference bandwidth, evaluated
% on 100 grid points from three bandwidths below the minimum to three above the
% maximum) and returns statistics summarizing its shape: the number of peaks, the
% height of the highest peak, the entropy, and two measures of asymmetry about
% the mean. Optionally, also counts the crossings of the curve through given
% heights, measures the area under the curve where it is lower than given
% heights, and the total variation of the curve within given distances of the
% mean. Heights and distances are in the units of the input, so the optional
% statistics depend on the scale of the data.
%
% ---INPUTS:
% x, the input data vector
% <can also produce additional outputs with the following optional settings>
% [opt] 'numcross': number of times the distribution crosses each given height
%           e.g., usage: DN_FitKernelSmooth(x,'numcross',[0.5,0.7]) for
%                        heights of 0.5 and 0.7
% [opt] 'area': area under the curve where it is below each given height.
%               Usage as for 'numcross' above
% [opt] 'arclength': total variation of the curve, sum(abs(diff(f)))*dx,
%               over the region within each given distance of the mean.
%               Usage as above.
%
% ---EXAMPLE USAGE:
% DN_FitKernelSmooth(x,'numcross',[0.05,0.1],'area',[0.1,0.2,0.4],'arclength',[0.5,1,2])
% returns all the basic outputs, plus those for numcross, area, and arclength
% for the thresholds given
%
% ---OUTPUTS:
% npeaks, the number of peaks (local maxima with a second difference below
%       -0.0002, i.e., clearly peaked)
% max, the height of the highest peak
% entropy, the entropy of the distribution, -sum(f*log(f)*dx), in nats
% asym, the probability mass above the mean divided by that below it (NaN if there
%       is essentially none below)
% plsym, the total variation of the curve below the mean divided by that
%       above the mean (NaN if there is none above)
% numcross_005, numcross_010, numcross_020, numcross_030, numcross_040,
% numcross_050, ...: the number of crossings of each threshold given to
%       'numcross' (named for the threshold to two decimal places, without the
%       point: 0.05 gives numcross_005)
% area_005, area_010, area_020, area_030, area_040, area_050, ...: the area
%       under the curve where it is below each threshold given to 'area'
% arclength_010, arclength_050, arclength_100, arclength_200, ...: the total
%       variation of the curve within each distance of the mean given to
%       'arclength', multiplied by the grid spacing

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

% -------------------------------------------------------------------------------
% Check inputs using the inputParser:
% -------------------------------------------------------------------------------
inputP = inputParser;
check_isPositive = @(x) validateattributes(x, {'numeric'}, {'positive'});
addParameter(inputP, 'area', [], check_isPositive);
addParameter(inputP, 'numcross', [], check_isPositive);
addParameter(inputP, 'arclength', [], check_isPositive);
parse(inputP, varargin{:});

% Make variables from results of input parser:
area = inputP.Results.area;
numcross = inputP.Results.numcross;
arclength = inputP.Results.arclength;
clear inputP;

% -------------------------------------------------------------------------------
% Preliminary definitions
m = mean(x);

% First compute the smoothed empirical distribution of values in the time series
[f, xi] = BF_KSDensity(x);
if any(isnan(f)) % constant data: no scale to smooth over
	out = NaN; return
end

% 1. Number of peaks
df = diff(f);
ddf = diff(df);
sdsp = ddf(BF_SignChange(df, 1));
out.npeaks = sum(sdsp < -0.0002); % 'large enough' maxima

% 2. Max
out.max = max(f); % maximum of the distribution

% 3. Entropy
out.entropy = -sum(f(f > 0) .* log(f(f > 0)) * (xi(2) - xi(1))); % entropy of the distribution

% 4. Assymetry
out1 = sum(f(xi > m) .* (xi(2) - xi(1)));
out2 = sum(f(xi < m) .* (xi(2) - xi(1)));
if out2 < 1e-10 % (essentially) no mass below the mean: the ratio is not meaningful
	out.asym = NaN;
else
	out.asym = out1 / out2;
end

% 5. Plsym
out1 = sum(abs(diff(f(xi < m))) .* (xi(2) - xi(1)));
out2 = sum(abs(diff(f(xi > m))) .* (xi(2) - xi(1)));
if out2 < 1e-10 % no variation above the mean
	out.plsym = NaN;
else
	out.plsym = out1 / out2;
end

% ------------------------------------------------------------------------------
% 6. Numcross
% ------------------------------------------------------------------------------
% Specified in input
if ~isempty(numcross) % calculate crossing statistics
	thresholds = numcross;
	for i = 1:length(thresholds)
		numCrosses = sum(BF_SignChange(f - thresholds(i)));
		outName = regexprep(sprintf('numcross_%.2f', thresholds(i)), '\.', ''); % remove dots from 2-d.pl.
		out.(outName) = numCrosses;
	end
end

% ------------------------------------------------------------------------------
% 7. Area
% ------------------------------------------------------------------------------
% Specified in input

if ~isempty(area) % calculate area statistics
	thresholds = area;
	for i = 1:length(thresholds)
		areaHere = sum(f(f < thresholds(i)) .* (xi(2) - xi(1))); % integral under this portion
		outName = regexprep(sprintf('area_%.2f', thresholds(i)), '\.', ''); % remove dots from 2-d.pl.
		out.(outName) = areaHere;
	end
end

% ------------------------------------------------------------------------------
% 8. Arc length
% ------------------------------------------------------------------------------
% Specified in input
if ~isempty(arclength) % calcualte arc length statistics
	thresholds = arclength;
	for i = 1:length(thresholds)
		% The integrand in the path length formula:
		fd = abs(diff(f(xi > m - thresholds(i) & xi < m + thresholds(i))));
		arclengthHere = sum(fd .* (xi(2) - xi(1)));
		outName = regexprep(sprintf('arclength_%.2f', thresholds(i)), '\.', ''); % remove dots from 2-d.pl.
		out.(outName) = arclengthHere;
	end
end

end
