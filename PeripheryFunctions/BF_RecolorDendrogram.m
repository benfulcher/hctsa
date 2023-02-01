function BF_RecolorDendrogram(h_dend);
% BF_RecolorDendrogram     Recolor a dendrogram using hctsa default colors
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

numLines = length(h_dend);

% Find (non-black) colored groups:
allColors = vertcat(h_dend.Color);
uniqueColors = unique(allColors,'rows');
isBlack = @(x) x(1)==0 && x(2)==0 && x(3)==0;
rowIsBlack = arrayfun(@(x)isBlack(uniqueColors(x,:)),1:size(uniqueColors,1));
uniqueColors(rowIsBlack,:) = [];

% Swap out weird colors for nicer hctsa default ones
numColors = size(uniqueColors,1);
colorCell = GiveMeColors(numColors);

for i = 1:numColors
    isMatch = zeros(numLines,1);
    for j = 1:numLines
        if allColors(j,:)==uniqueColors(i,:)
            isMatch(j) = 1;
        end
    end
    matchInd = find(isMatch);
    numMatches = length(matchInd);
    for j = 1:numMatches
        h_dend(matchInd(j)).Color = colorCell{i};
    end
end

end
