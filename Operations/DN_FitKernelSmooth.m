% DN_FitKernelSmooth
% 
% Fits a kernel-smoothed distribution to the data using the ksdensity function
% from Matlab's Statistics Toolbox and returns a set of simple statistics.
% 
% INPUTS:
% x, the input time series
% <can also produce additional outputs with the following optional settings>
% [opt] 'numcross': number of times the distribution crosses the given threshold
%           e.g., usage: DN_FitKernelSmooth(x,'numcross',[0.5,0.7]) for
%                        thresholds of 0.5 and 0.7
% [opt] 'area': area under where the distribution crosses the given thresholds.
%               Usage as for 'numcross' above
% [opt] 'arclength': arclength between where the distribution passes given
%       thresholds. Usage as above.
% 
% EXAMPLE USAGE:                  
% DN_FitKernelSmooth(x,'numcross',[0.05,0.1],'area',[0.1,0.2,0.4],'arclength',[0.5,1,2])
% returns all the basic outputs, plus those for numcross, area, and arclength
% for the thresholds given
% 
% Outputs are set of statistics summarizing the obtained distribution, including
% the number of peaks, the distributional entropy, the number of times the curve
% crosses fifixed probability thresholds, the area under the curve for fifixed
% probability thresholds, the arc length, and the symmetry of probability
% density above and below the mean.
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

function out = DN_FitKernelSmooth(x,varargin)
% Ben Fulcher, 2009

% Preliminary definitions
m = mean(x);

% First compute the smoothed empirical distribution of values in the time series
[f, xi] = ksdensity(x);

% 1. Number of peaks
df = diff(f);
ddf = diff(df);
sdsp = ddf(BF_sgnchange(df,1));
out.npeaks = sum(sdsp < -0.0002); % 'large enough' maxima

% 2. Max
out.max = max(f); % maximum of the distribution

% 3. Entropy
out.entropy = - sum(f(f > 0).*log(f(f > 0))*(xi(2)-xi(1))); % entropy of the distribution

% 4. Assymetry
out1 = sum(f(xi > m).*(xi(2)-xi(1)));
out2 = sum(f(xi < m).*(xi(2)-xi(1)));
out.asym = out1/out2;

% 5. Plsym
out1 = sum(abs(diff(f(xi < m))).*(xi(2)-xi(1)));
out2 = sum(abs(diff(f(xi > m))).*(xi(2)-xi(1)));
out.plsym = out1/out2;

% 6. Numcross
% Specified in input
varhere = strcmp(varargin,'numcross');
if any(varhere) % calculate crossing statistics
    thresholds = varargin{find(varhere,1,'first') + 1};
    for i = 1:length(thresholds)
        ncrosses = sum(BF_sgnchange(f - thresholds(i)));
        outname = regexprep(sprintf('numcross_%.2f',thresholds(i)),'\.',''); % remove dots from 2-d.pl.
        eval(sprintf('out.%s = ncrosses;',outname));
    end
end

% 7. Area
% Specified in input
varhere = strcmp(varargin,'area');
if any(varhere) % calculate area statistics
    thresholds = varargin{find(varhere,1,'first') + 1};
    for i = 1:length(thresholds)
        areahere = sum(f(f < thresholds(i)).*(xi(2)-xi(1))); % integral under this portion
        outname = regexprep(sprintf('area_%.2f',thresholds(i)),'\.',''); % remove dots from 2-d.pl.
        eval(sprintf('out.%s = areahere;',outname));
    end
end

% 8. Arc length
% Specified in input
varhere = strcmp(varargin,'arclength');
if any(varhere) % calcualte arc length statistics
    thresholds = varargin{find(varhere,1,'first') + 1};
    for i = 1:length(thresholds)
        % The integrand in the path length formula:
        fd = abs(diff(f(xi > m - thresholds(i) & xi < m + thresholds(i))));
        arclengthhere = sum(fd.*(xi(2)-xi(1)));
        outname = regexprep(sprintf('arclength_%.2f',thresholds(i)),'\.',''); % remove dots from 2-d.pl.
        eval(sprintf('out.%s = arclengthhere;',outname));
    end
end


end