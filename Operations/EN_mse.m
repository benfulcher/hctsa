function out = EN_mse(y,scaleRange,m,r,preProcessHow)
% EN_mse multiscale entropy for a time series
%
% As per "Multiscale entropy analysis of biological signals",
% Costa, Goldberger and Peng, PRE, 71, 021906 (2005)
% http://physionet.comp.nus.edu.sg/physiotools/mse/papers/pre-2005.pdf
%
%---INPUTS:
% scaleRange: a vector of scales (default: 1:10)
% m: embedding dimension/length of sequence to match (default: 2)
% r: similarity threshold for matching (default: 0.15)
% preProcessHow: how to preprocess the data (default: do not)
%
%
% Original C implementation and docs here:
% http://physionet.org/physiotools/mse/tutorial/node3.html

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
% Check inputs, set defaults
%-------------------------------------------------------------------------------
if nargin < 2
    m = 2;
end
if nargin < 3
    r = 0.15;
end
if nargin < 4
    scaleRange = 1:10;
end
if nargin < 5
    preProcessHow = '';
end

%-------------------------------------------------------------------------------
% Impose a minimum time-series length of 20 samples to perform a SampEn
% (should probably be even higher...?)
minTSLength = 20;

doPlot = 0; % whether to plot outputs
numScales = length(scaleRange);

%-------------------------------------------------------------------------------
% Preprocess
%-------------------------------------------------------------------------------
% Do the specified pre-processing BEFORE applying the coarse-graining
if ~isempty(preProcessHow)
    y = zscore(BF_preprocess(y,preProcessHow));
end

%-------------------------------------------------------------------------------
% Coarse-graining across scales:
%-------------------------------------------------------------------------------
% cf. Eq. (16) in Costa et al. (2005)
y_cg = cell(numScales,1);
for i = 1:numScales
    % Want non-overlapping windows of length scaleRange(i)
    bufferSize = scaleRange(i);
    y_buffer = BF_makeBuffer(y,bufferSize);
    y_cg{i} = mean(y_buffer,2);
end

%-------------------------------------------------------------------------------
% Run sample entropy for each m and r value at each scale:
%-------------------------------------------------------------------------------
sampEns = zeros(numScales,1);
for si = 1:numScales
    if length(y_cg{si}) >= minTSLength
        sampEnStruct = EN_SampEn(y_cg{si},m,r);
        sampEns(si) = sampEnStruct.(sprintf('sampen%u',m));
    else
        sampEns(si) = NaN;
    end
end

%-------------------------------------------------------------------------------
% Outputs: multiscale entropy
%-------------------------------------------------------------------------------
if all(isnan(sampEns))
    if ~isempty(preProcessHow)
        ppText = sprintf('after %s pre-processing',preProcessHow);
    else
        ppText = '';
    end
    warning('Not enough samples (%u %s) to compute SampEn at multiple scales',...
                    length(y),ppText)
    out = NaN; return
end

if doPlot
    figure('color','w')
    subplot(2,1,1);
    plot(y);
    subplot(2,1,2);
    plot(sampEns,'o-k')
end

% Output raw values
for i = 1:numScales
    out.(sprintf('sampen_s%u',scaleRange(i))) = sampEns(i);
end

%-------------------------------------------------------------------------------
% Summary statistics of the variation:
%-------------------------------------------------------------------------------
% Maximum, and where it occurred
[out.maxSampEn,maxInd] = nanmax(sampEns);
out.maxScale = scaleRange(maxInd);
% Minimum, and where it occurred
[out.minSampEn,minInd] = nanmin(sampEns);
out.minScale = scaleRange(minInd);
% Mean, std, coefficient of variation:
out.meanSampEn = nanmean(sampEns);
out.stdSampEn = nanstd(sampEns);
out.cvSampEn = out.stdSampEn/out.meanSampEn;
% Mean change across the range of scales:
out.meanch = nanmean(diff(sampEns));

end
