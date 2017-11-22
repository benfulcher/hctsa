function out = SC_MMA(y,doOverlap,scaleRange,qRange)
% SC_MMA   Physionet implementation of multiscale multifractal analysis
%
% Scale-dependent estimates of multifractal scaling in a time series.

% ------------------------------------------------------------------------------
% Modified by Ben Fulcher for use in hctsa, 2015-05-12.
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
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Copyright (C) 2014 Jan Gieraltowski
% ------------------------------------------------------------------------------
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
% PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 59 Temple
% Place - Suite 330, Boston, MA 02111-1307, USA.
%
% Author: Jan Gieraltowski
% Warsaw University of Technology
% Faculty of Physics
% gieraltowski@if.pw.edu.pl
% http://gieraltowski.fizyka.pw.edu.pl/
%
% Method was first proposed in:
% J. Gieraltowski, J. J. Zebrowski, and R. Baranowski,
% Multiscale multifractal analysis of heart rate variability recordings
% with a large number of occurrences of arrhythmia,
% Phys. Rev. E 85, 021915 (2012).
% http://dx.doi.org/10.1103/PhysRevE.85.021915
%
% Please cite the above publication when referencing this material,
% and also include the standard citation for PhysioNet:
% Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG,
% Mietus JE, Moody GB, Peng C-K, Stanley HE.
% PhysioBank, PhysioToolkit, and PhysioNet: Components of a New Research Resource
% for Complex Physiologic Signals.
% Circulation 101(23):e215-e220
% [Circulation Electronic Pages; http://circ.ahajournals.org/cgi/content/full/101/23/e215]; 2000 (June 13).
% ------------------------------------------------------------------------------

% Generate a plot:
doPlot = false;

% Time-series length:
N = length(y);

% ------------------------------------------------------------------------------
% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(doOverlap)
    doOverlap = false;
    % 0 - time series is partitioned into non overlapping windows of analysis,
    % 1- time series is partitioned into overlapping windows of analysis, step between consecutive windows is = 1 (much longer calculations)
end

if nargin < 3 || isempty(scaleRange)
    scaleRange = [10,round(N/40)];
    % minimal s scale used, when calculating Fq(s) functions family (default 10)
    % maximal s scale used, when calculating Fq(s) functions family, has to be multiple of 5 (default 600; in general should be near to N/50, where N is a time series length)
end
minScale = scaleRange(1);
maxScale = scaleRange(2);
if (maxScale/5) < minScale
    warning('Time-series (N=%u) too short for multiscale multifractal analysis',N);
    out = NaN;
    return
elseif rem(maxScale,5)~=0
    maxScale = round(maxScale/5)*5;
    fprintf(1,'adjusted maxScale to %u\n',maxScale);
end

if nargin < 4 || isempty(qRange)
    qRange = [-5,5];
    % minimal/maximal multifractal parameter q used (default -5)
end
qMin = qRange(1);
qMax = qRange(2);

qList = qMin:0.1:qMax;
qList(qList == 0) = 0.0001;

% if nargin < 5
%     precisionMode = 0;
%     % 0 (default) better looking plot, smaller files, faster calculations;
%     % set 1 for enhanced precision (smaller q and s steps)
% end

% ------------------------------------------------------------------------------

signal = y;

prof = cumsum(signal);
slength = size(prof);
fqs = [];

numIncrements = 20;
sListFull = unique(round(linspace(minScale,maxScale,numIncrements)));

timer = tic;
for s = sListFull

    if doOverlap
        vec = [0:s-1];
        ind = [1:slength-s+1]';
        coordinates = bsxfun(@plus, vec, ind);
    else
        ind2 = [1:size(prof,1)];
        coordinates = reshape(ind2(1:(size(prof,1)-mod(size(prof,1),s))),s,(size(prof,1)-mod(size(prof,1),s))/s)';
    end

    segments = prof(coordinates);
    xbase = [1:1:s];
    f2nis = [];

    for ni = 1:size(segments,1)
        seg = segments(ni,:);
        fit = polyfit(xbase,seg,2);
        variance = mean((seg - polyval(fit,xbase)).^2);
        f2nis(end+1) = variance;
    end

    for q = qList
        fqs = [fqs; q s (mean(f2nis.^(q/2)))^(1/q)];
    end
end
% fprintf(1,'Detrended fluctuations computed in %s\n',BF_thetime(toc(timer)));

fqsll = [fqs(:,1) fqs(:,2) log(fqs(:,2)) log(fqs(:,3))];


% ------------------------------------------------------------------------------
% Now compute hurst exponents as the gradients of F(q) curves
% ------------------------------------------------------------------------------
hqs = [];

% if precisionMode == 0
%     % Take 11 points through the space
%     sspacing = ((maxScale/5)-minScale)/10;
%     sList = minScale:sspacing:(maxScale/5);
%     % Coarser sampling of q space
%     qList = qMin:1:qMax; qList(qList == 0) = 0.0001;
% else
    % sspacing = 1;
% sList = minScale:sspacing:(maxScale/5);

if sum(sListFull<=maxScale/5)>=10
    sList = sListFull(sListFull<=maxScale/5);
else
    % sample higher in the scale dimension:
    sspacing = ((maxScale/5)-minScale)/10;
    sList = minScale:sspacing:(maxScale/5);
end

% Coarser sampling of q space
qList = qMin:0.5:qMax; qList(qList == 0) = 0.0001;
% end

% sList = minScale:sspacing:(maxScale/5);

hqs = zeros(length(qList),length(sList));
for si = 1:length(sList)
    sit = sList(si);
    for qi = 1:length(qList)
        qit = qList(qi);

        fitTemp = fqsll(fqsll(:,1) == qit & fqsll(:,2) >= sit & fqsll(:,2) <= 5*sit,:);
        hTemp = polyfit(fitTemp(:,3),fitTemp(:,4),1);

        hqs(qi,si) = hTemp(1);
        % hqs = [hqs; qit 3*sit hTemp(1)];
    end
end

sListScaled = sList*3; % Not completely on top of the algorithm, but for some reason
                       % this was recorded as a multiple of 3 in the original algorithm

% hqsPlotData = reshape(hqs(:,3),size(qList,2),length(minScale:sspacing:(maxScale/5)));

% ------------------------------------------------------------------------------
% Plotting
% ------------------------------------------------------------------------------

if doPlot
    if max(max(hqs)) < 1.5
        hLim = 1.5;
    elseif max(max(hqs)) < 2.5
        hLim = 2.5;
    else
        hLim = ceil((max(max(hqs))*10))/10;
    end

    f = figure('color','w'); box('on');
    hqsplot = surf(sListScaled,qList,hqs);
    colormap(jet);
    colorbar;
    caxis([0,hLim]);
    set(gca,'YDir','reverse');
    view(-62,50);
    axis([sListScaled(1),sListScaled(end),qMin,qMax,0,hLim]);
    xlabel('scale')
    ylabel('q')
    zlabel('h')
end

% ------------------------------------------------------------------------------
% Output statistics:
% ------------------------------------------------------------------------------

% hqsPlotData (scale,q)

% Global properties:
allExponents = hqs(:);
out.meanHurstExponent = mean(allExponents);
out.stdHurstExponent = std(allExponents);
out.minHurstExponent = min(allExponents);
out.maxHurstExponent = max(allExponents);

% Changes with scale:
out.scaleHurstStd = std(mean(hqs,1));
out.scaleHurstTrend = GiveMeGradient(sListScaled,mean(hqs,1));

% Changes with q:
out.qHurstStd = std(mean(hqs,2));
out.qHurstTrend = GiveMeGradient(qList,mean(hqs,2));

% max/min points are where in scale/q space?
[qi,si] = find(hqs==max(hqs(:)),1);
out.maxHurstQ = qList(qi);
out.maxHurstScale = sListScaled(si);
[qi,si] = find(hqs==min(hqs(:)),1);
out.minHurstQ = qList(qi);
out.minHurstScale = sListScaled(si);

% phase transitions
% there is some peak or trough somewhere, so the standard deviation across
% scales or q is inconsistent
stdS = std(hqs,[],1);
stdQ = std(hqs,[],2);
out.stdStdHurstQ = std(stdQ); % will be large if variance changes alot with Q
out.stdStdHurstScale = std(stdS); % will be large if variance changes alot with scale

% saveas(hqsplot, [dirname filesep 'MMA_results' filesep 'MMA_' n1 '.jpg']);

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
