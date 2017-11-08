function out = NL_DVV(x,m,numDVs,nd,Ntv,numSurr,randomSeed)
% NL_DVV 	Delay Vector Variance method for real and complex signals.
%
% Uses predictability of the signal in phase space to characterize the
% time series.
%
% This function uses the original code from the DVV toolbox to do the computation
% and produces statistics on the outputs -- comparing the DVV curves for both
% the real and the surrogate data.
%
%---USAGE:
% outputStats = NL_DVV(x, m, numDVs, nd, Ntv, numSurr, randomSeed)
%
%---INPUTS:
% x:            Original real-valued or complex time series
% m:            Delay embedding dimension
% Ntv:          Number of points on horizontal axes
% numDVs:	    Number of reference DVs to consider
% nd:           Span over which to perform DVV
% Ntv:          Number of points on the horizontal axis
% numSurr:      Number of surrogates to compare to
% randomSeed:   How to control the random seed for reproducibility
%
% A Delay Vector Variance (DVV) toolbox for MATLAB
% (c) Copyright Danilo P. Mandic 2008
% http://www.commsp.ee.ic.ac.uk/~mandic/dvv.htm
% http://www.commsp.ee.ic.ac.uk/~mandic/dvv/papers/dvv_proj.pdf
%
% ------------------------------------------------------------------------------
% Modified by Ben Fulcher, 2015-05-13, for use in hctsa.
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
% ------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the Free
%   Software Foundation; either version 2 of the License, or (at your option)
%   any later version.
%
%   This program is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
%   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
%   more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to Free Software
%   Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot output plot:
doPlot = false;

% Talk to me:
beVocal = false;

% ------------------------------------------------------------------------------
% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 1
	error('Not enough input arguments');
end
if nargin < 2 || isempty(m)
	m = 3;
end
if nargin < 3 || isempty(numDVs)
	numDVs = 100;
end
if nargin < 4 || isempty(nd)
	nd = 2.0;
end
if nargin < 5 || isempty(Ntv)
	Ntv = 25*nd;
end
if nargin < 6 || isempty(numSurr)
    numSurr = 10;
end
if nargin < 7
    randomSeed = [];
end

% ------------------------------------------------------------------------------
% Preliminaries
% ------------------------------------------------------------------------------
BF_ResetSeed(randomSeed); % Reset the random seed if specified

% ------------------------------------------------------------------------------
% Run DVV on input data:
% ------------------------------------------------------------------------------
if beVocal, fprintf(1,'dvv on data...'); end
dvv_data = DVV_dvv(x,m,numDVs,nd,Ntv);
if beVocal, fprintf(1,' Done.\n'); end

% ------------------------------------------------------------------------------
% Generate surrogate data and run DVV on the surrogates
% ------------------------------------------------------------------------------
if beVocal, fprintf(1,'Generating %u surrogates...',numSurr); end
x_surr = DVV_surrogate(x, numSurr);
if beVocal, fprintf(1,' Done.\n'); end

if beVocal, fprintf(1,'Computing dvv for surrogates...'); end
dvv_surr = zeros(Ntv,2,numSurr);
for i = 1:numSurr
    dvv_surr(:,:,i) = DVV_dvv(x_surr(:,i),m,numDVs,nd,Ntv);
end
mean_dvv_surr = mean(dvv_surr,3);
if beVocal, fprintf(1,' Done.\n'); end

% ------------------------------------------------------------------------------
% Plotting:
% ------------------------------------------------------------------------------
if doPlot
    f = figure('color','w'); box('on');
    subplot(1,2,1); hold on
    plot(dvv_data(:,1),dvv_data(:,2),'.-k')
    plot(mean_dvv_surr(:,1),mean_dvv_surr(:,2),'.-r')
    legend({'data','surrogate mean'})
    xlabel('~distance')
    ylabel('target variance')

    % DVV scatter plot
    subplot(1,2,2); hold on
    biSector = 0:0.5:1;
    plot(biSector, biSector, ':k'); axis([ 0 1 0 1]);
    title('DVV Scatter Plot'); xlabel('Original'); ylabel('Surrogates'); grid on; hold on
    errorbar(dvv_data(:,2),mean_dvv_surr(:,2),std(dvv_surr(:,2,:),[],3));
end


% ------------------------------------------------------------------------------
% Quantify and output patterns in the distance versus target variance data
% (for both real data and surrogates)
% ------------------------------------------------------------------------------

dvv_data_notNaN = dvv_data(~isnan(dvv_data(:,2)),:);

% Is there a trend in the data?
out.trend = GiveMeGradient(dvv_data_notNaN(:,1),dvv_data_notNaN(:,2));

% What's the distribution of values like
out.max = max(dvv_data_notNaN(:,2));
out.min = min(dvv_data_notNaN(:,2));
out.mean = mean(dvv_data_notNaN(:,2));

% Analysis of successive differences in the real data:
% (assume that NaNs occur in contiguous block at the start; for low distances)
out.meanDiff = mean(diff(dvv_data_notNaN(:,2)));
out.stdDiff = std(diff(dvv_data_notNaN(:,2)));
out.trendDiff = GiveMeGradient(dvv_data_notNaN(1:end-1,1),diff(dvv_data_notNaN(:,2)));

% Are there large differences between real and surrogate data?
notNaN = ~isnan(dvv_data(:,2)) & ~isnan(mean_dvv_surr(:,2));
out.rmsDiffSurr = sqrt(mean((dvv_data(notNaN,2) - mean_dvv_surr(notNaN,2)).^2));
out.meanDiffSurr = mean(dvv_data(notNaN,2) - mean_dvv_surr(notNaN,2));

% What's the correlation between the original and surrogate DVV values
out.dataSurrCorr = corr(dvv_data(notNaN,2),mean_dvv_surr(notNaN,2));

% Trend in DVV scatter plot
out.trendDataSurr = GiveMeGradient(dvv_data(notNaN,2),mean_dvv_surr(notNaN,2));

% Number of zero crossings of data/surrogate data curve
out.numZeroCrossings = sum(BF_sgnchange(dvv_data(notNaN,2)-mean_dvv_surr(notNaN,2)));

% Difference in trend between real and surrogate:
out.trendSurr = GiveMeGradient(mean_dvv_surr(notNaN,1),mean_dvv_surr(notNaN,2));
out.meanDiffTrendSurr = out.trendSurr - out.trend;

% Mean probability assuming Gaussian distribution of test stat for surrogates:
probs = zeros(sum(notNaN),1);
fnotNaN = find(notNaN);
for i = 1:sum(notNaN)
    probs(i) = normcdf(dvv_data(fnotNaN(i),2),mean(dvv_surr(fnotNaN(i),2,:)),std(dvv_surr(fnotNaN(i),2,:)));
end
out.meanNormCDF = mean(probs);


% ------------------------------------------------------------------------------
function m = GiveMeGradient(xData,yData)
    if size(xData,1) ~= size(yData,1);
        yData = yData';
    end
    p = polyfit(xData,yData,1);
    m = p(1);
end
% ------------------------------------------------------------------------------

end
