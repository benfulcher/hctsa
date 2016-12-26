function out = EN_Randomize(y,randomizeHow,randomSeed)
% EN_Randomize  How time-series properties change with increasing randomization.
%
% Progressively randomizes the input time series according to a specified
% randomization procedure
%
% The procedure is repeated 2N times, where N is the length of the time series.
%
%---INPUTS:
% y, the input (z-scored) time series
%
% randomizeHow, specifies the randomization scheme for each iteration:
%      (i) 'statdist' -- substitutes a random element of the time series with
%                           one from the original time-series distribution
%      (ii) 'dyndist' -- overwrites a random element of the time
%                       series with another random element
%      (iii) 'permute' -- permutes pairs of elements of the time
%                       series randomly
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed
%
%---OUTPUTS: summarize how the properties change as one of these
% randomization procedures is iterated, including the cross correlation with the
% original time series, the autocorrelation of the randomized time series, its
% entropy, and stationarity.
%
% These statistics are calculated every N/10 iterations, and thus 20 times
% throughout the process in total.
%
% Most statistics measure how these properties decay with randomization, by
% fitting a function f(x) = Aexp(Bx).

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
% Check toolboxes, and a z-scored time series:
% ------------------------------------------------------------------------------

% Check a curve-fitting toolbox license is available:
BF_CheckToolbox('curve_fitting_toolbox');

% Check input time series is z-scored:
if ~BF_iszscored(y)
    warning('The input time series should be z-scored for EN_Randomize.')
end

% ------------------------------------------------------------------------------
%% Check inputs:
% ------------------------------------------------------------------------------
% randomizeHow, how to do the randomization:
if nargin < 2 || isempty(randomizeHow)
    randomizeHow = 'statdist'; % use statdist by default
end

% randomSeed: how to treat the randomization
if nargin < 3
    randomSeed = []; % default
end
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Preliminaries
% ------------------------------------------------------------------------------

doPlot = 0; % Don't plot to screen by default:
N = length(y); % length of the time series

% Set up the points through the randomization process at which the
% calculation of stats will occur:
randp_max = 2; % time series has been randomized to double its length
rand_inc = 0.1; % this proportion of the time series has been randomized between calculations
numCalcs = randp_max/rand_inc; % number of calculations required
calc_ints = floor(randp_max*N/numCalcs);
calc_pts = (0:calc_ints:randp_max*N);
if calc_pts(end) ~= randp_max*N;
    calc_pts = [calc_pts,randp_max*N];
end
numCalcs = length(calc_pts); % some rounding issues inevitable

statNames = {'xcn1', 'xc1', 'd1', 'ac1', 'ac2', 'ac3', 'ac4', 'sampen2_015', 'statav5', 'swss5_1'};
numStats = length(statNames);
stats = zeros(numCalcs,numStats); % record a stat at each randomization increment

y_rand = y; % this vector will be randomized

stats(1,:) = CalculateStats(y,y_rand); % initial condition: apply on itself

% Control the random seed (for reproducibility):
BF_ResetSeed(randomSeed);

%-------------------------------------------------------------------------------
% Do the randomization
%-------------------------------------------------------------------------------
% fprintf(1,'%u/%u calculation points',numCalcs,N*randp_max)

for i = 1:N*randp_max
    switch randomizeHow

        case 'statdist'
            % randomize by substituting a random element of the time series by
            % a random element from the static original time series distribution
            y_rand(randi(N)) = y(randi(N));

        case 'dyndist'
            % randomize by substituting a random element of the time series
            % by a random element of the current, already partially randomized,
            % time series
            y_rand(randi(N)) = y_rand(randi(N));

        case 'permute'
            % randomize by swapping elements of the time series so that
            % the distribution remains static; only temporal properties will change
            randis = randi(N,[2, 1]);
            y_rand(randis(1)) = y_rand(randis(2));
            y_rand(randis(2)) = y_rand(randis(1));

        otherwise
            error('Unknown randomization method ''%s''.',randomizeHow);
    end

    if any(calc_pts==i)
        stats(calc_pts == i,:) = CalculateStats(y,y_rand);
    end

end
% fprintf(1,'Randomization took %s',BF_thetime(toc(randTimer)));

if doPlot
    f = figure('color','w'); box('on');
    plot(stats,'.-');
end

% ------------------------------------------------------------------------------
%% Fit exponentials to outputs:
% ------------------------------------------------------------------------------
r = (1:size(stats,1))'; % gives an 'x-axis' for fitting

% 1) xcn1: cross correlation at lag of -1
% 2) xc1: cross correlation at lag 1
% 3) d1: norm of differences between original and randomized time series
% 4) ac1
% 5) ac2
% 6) ac3
% 7) ac4
% 8) sample entropy
% 9) statav5
% 10) swss5_1

out = struct();
for i = 1:length(statNames)
    % Exponential fits:
    switch statNames{i}
    case {'xcn1','xc1'}
        startPoint = [stats(1,i),-0.1];
        [c,gof] = f_fix_exp(r,stats(:,i),startPoint,0);
    case {'ac1','ac2','ac3'}
        startPoint = [stats(1,i),-0.2];
        [c,gof] = f_fix_exp(r,stats(:,i),startPoint,0);
    case 'ac4'
        startPoint = [stats(1,i),-0.4];
        [c,gof] = f_fix_exp(r,stats(:,i),startPoint,0);
    case {'d1','sampen2_015'}
        startPoint = [-stats(end,i),-0.2,stats(end,i)];
        [c,gof] = f_fix_exp(r,stats(:,i),startPoint,1);
    case {'statav5','swss5_1'}
        startPoint = [-stats(end,i),-0.1,stats(end,i)];
        [c,gof] = f_fix_exp(r,stats(:,i),startPoint,1);
    end
    out = assignExpStats(out,c,gof,statNames{i});

    % Extra statistics:
    out = assignExtraStats(out,stats(:,i),statNames{i});
end

% ------------------------------------------------------------------------------
    function out = CalculateStats(y,y_rand)
        % Calculate statistics comparing a time series, y, and a randomized
        % version of it, y_rand

        % Cross Correlation to original signal
        xc = xcorr(y,y_rand,1,'coeff');
        xcn1 = xc(1);
        xc1 = xc(3);

        % Norm of differences between original and randomized signals
        d1 = norm(y - y_rand) / length(y);

        % Autocorrelation
        autoCorrs = CO_AutoCorr(y_rand,1:4,'Fourier');
        ac1 = autoCorrs(1);
        ac2 = autoCorrs(2);
        ac3 = autoCorrs(3);
        ac4 = autoCorrs(4);

        % 2-bit LZ complexity:
        % LZcomplex = EN_MS_LZcomplexity(y,3);

        % SampEn(2,0.2,1):
        sampenStruct = EN_SampEn(y_rand,2,0.15);
        sampen2_015 = sampenStruct.quadSampEn2;

        % Stationarity
        statav5 = SY_StatAv(y_rand,'seg',5);
        swss5_1 = SY_SlidingWindow(y_rand,'std','std',5,1);

        out = [xcn1, xc1, d1, ac1, ac2, ac3, ac4, sampen2_015, statav5, swss5_1];
    end
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
function thehp = SUB_gethp(v)
	if v(end) > v(1)
		thehp = find(v > 0.5*(v(end)+v(1)),1,'first');
	else
		thehp = find(v < 0.5*(v(end)+v(1)),1,'first'); % last?
	end
end
% ------------------------------------------------------------------------------
function [c,gof] = f_fix_exp(r,dataVector,startPoint,addOffset)
    % Fits an exponential to the data vector across data points r

    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',startPoint);
    if addOffset
        f = fittype('a*exp(b*x)+c','options',s);
        f_x = @(c,x) c.a*exp(c.b*x)+c.c;
    else
        f = fittype('a*exp(b*x)','options',s);
        f_x = @(c,x) c.a*exp(c.b*x);
    end
    try
        [c, gof] = fit(r,dataVector,f);
    catch
        warning('Exponential fit failed :(')
        c = struct('a',NaN,'b',NaN,'c',NaN);
        gof = struct('rsquare',NaN,'rmse',NaN);
    end
    if doPlot;
        figure('color','w'); hold on;
        plot(r,dataVector,'x-k');
        xr = linspace(min(r),max(r),100);
        plot(xr,f_x(c,xr))
    end
end
%-------------------------------------------------------------------------------
function out = assignExpStats(out,c,gof,fieldName)
    % Assigns relevant stats from an exponential fit result, [c,gof]
    out.([fieldName,'fexpa']) = c.a;
    out.([fieldName,'fexpb']) = c.b;
    if (isstruct(c) && isfield(c,'c')) || ismember('c',coeffnames(c))
        out.([fieldName,'fexpc']) = c.c;
    end
    out.([fieldName,'fexpr2']) = gof.rsquare;
    out.([fieldName,'fexprmse']) = gof.rmse;
end
%-------------------------------------------------------------------------------
function out = assignExtraStats(out,dataVector,fieldName)
    % Assigns 2 extra statistics about a data vector:
    out.([fieldName,'diff']) = abs((dataVector(end)-dataVector(1))/dataVector(1));
    out.([fieldName,'hp']) = SUB_gethp(dataVector);
end

end
