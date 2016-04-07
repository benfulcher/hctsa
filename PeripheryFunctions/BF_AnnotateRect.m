function BF_AnnotateRect(whatCfn,featureVector,groupLabels,numClasses,colors,ax)
% BF_AnnotateRect Annotates rectangles under distributions for group classification
%
%---INPUTS:
% whatCfn, the classification method to use
% featureVector, the values given to each observation
% groupLabels, the class labels assigned to each observation
% numClasses, the number of classes (usually the max of the groupLabels)
% colors, a custom cell of colors for plotting
% ax, the handle for the axes to annotate

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

if isempty(whatCfn)
    whatCfn = 'fast_linear';
end
if nargin < 4
    numClasses = max(groupLabels);
end
if nargin < 5
    colors = BF_getcmap('dark2',numClasses,1);
end
if nargin < 6
    ax = gca;
end

% How many increments across the axis to determine classification labels:
numIncrements = 100;

%--------------------------------------------------------------------
% Learn a basic classification model from the 1-d feature vector:
[~,Mdl] = GiveMeCfn(whatCfn,featureVector,groupLabels,...
                    featureVector,groupLabels,numClasses,0,'balancedAcc');
xPlotted = ax.XLim;
predRange = linspace(xPlotted(1),xPlotted(2),numIncrements)'; % numIncrements-grid through x
predLabels = predict(Mdl,predRange);

%-------------------------------------------------------------------------------
% Annotate rectangles under the distribution reflecting the predictive model:
yPlotted = ax.YLim;
rectHeight = 0.1*diff(ax.YLim);
ind = 1;
while ind < numIncrements
    rectLength = find(predLabels(ind+1:end)~=predLabels(ind),1,'first');
    if isempty(rectLength)
        rectLength = length(predLabels)-ind;
    end
    rectangle('Position',[predRange(ind),-rectHeight*1.1,predRange(ind+rectLength)-predRange(ind),rectHeight],...
                                'FaceColor',colors{predLabels(ind)});
    ind = ind+rectLength;
end
ax.YLim(1) = -rectHeight*1.1;
ax.XLim = xPlotted;
ax.YLim(2) = yPlotted(2);

end
