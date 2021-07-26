function [N,binEdges] = BF_SimpleBinner(xData,numBins)
% BF_SimpleBinner Generate a histogram from equally spaced bins.
% An alternative to the more opaque method of Matlab's internal histcounts function.
%
%---INPUTS:
% xData, a data vector.
% numBins, the number of bins.
%---OUTPUTS:
% N, the counts
% binEdges, the extremeties of the bins.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

minX = min(xData);
maxX = max(xData);

% Linearly spaced bins:
binEdges = linspace(minX,maxX,numBins+1);

N = zeros(numBins,1);
for i = 1:numBins
    if i < numBins
        N(i) = sum(xData >= binEdges(i) & xData < binEdges(i+1));
    else
        % the final bin
        N(i) = sum(xData >= binEdges(i) & xData <= binEdges(i+1));
    end
end

end
