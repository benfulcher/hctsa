%%
%%  Java Information Dynamics Toolkit (JIDT)
%%  Copyright (C) 2012, Joseph T. Lizier
%%  
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%  
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%  
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%

% function h = plotLocalInfoValues(localResults, plotOptions)
%
% J. Lizier, 2012.
%
% Plots the local values, both the positive and negative, in blues and reds respectively, on the current figure.
%
% Inputs:
% - localResults - local values to be plotted. Can be native octave or java array
% - plotOptions - structure (optional) containing the following variables:
%   - plotRows - how many rows to plot (default is all)
%   - plotCols - how many columns to plot (default is all)
%   - plotStartRow - which row to start plotting from (default is 1)
%   - plotStartCol - which column to start plotting from (default is 1)
%   - scaleColoursToSubsetOfPlot - whether to plot the darkest blue and red for the max and min of
%      the entire profile, or use the darkest for the extremes within the subset that we are plotting (default true)
%   - scaleColoursToExtremes - stretch the darkest red and blue to the max and min of the local values, 
%      regardless of how imbalanced the max and mins are. (Default is false.)
%   - mainSignVectorLength - length of the longest component (positive or negative values) for the colourmap (default 1024)
%   - scalingMainComponent - what proportion to make the darkest shade of the primary colour (default .25)
%   - scalingScdryComponent - what proportion to make the lighter shades with green added (default .4)
%   - gammaPower - scale the colours non-linearly with this exponent (default 1 - linear scaling)
%
% Outputs:
% - h - return value of imagesc
%
function h = plotLocalInfoValues(localResults, plotOptions)

	% Set the colormap to have blue for positive, red for negative, scaled to the max and min of our local values

	if (nargin < 2)
		plotOptions = {};
	end
	if not(isfield(plotOptions, 'plotRows'))
		plotOptions.plotRows = size(localResults, 1);
	end
	if not(isfield(plotOptions, 'plotCols'))
		plotOptions.plotCols = size(localResults, 2);
	end
	if not(isfield(plotOptions, 'plotStartRow'))
		plotOptions.plotStartRow = 1;
	end
	if not(isfield(plotOptions, 'plotStartCol'))
		plotOptions.plotStartCol = 1;
	end
	if not(isfield(plotOptions, 'scaleColoursToSubsetOfPlot'))
		plotOptions.scaleColoursToSubsetOfPlot = true;
	end
	if not(isfield(plotOptions, 'scaleColoursToExtremes'))
		plotOptions.scaleColoursToExtremes = false;
	end
	if not(isfield(plotOptions, 'mainSignVectorLength'))
		plotOptions.mainSignVectorLength = 1024;
	end
	if not(isfield(plotOptions, 'scalingMainComponent'))
		plotOptions.scalingMainComponent = .25;
	end
	if not(isfield(plotOptions, 'scalingScdryComponent'))
		plotOptions.scalingScdryComponent = .4;
	end
	if not(isfield(plotOptions, 'gammaPower'))
		plotOptions.gammaPower = 1;
	end
	if (plotOptions.plotRows > 1000)
		fprintf('*** Limiting number of plotted rows to 1000\n');
		plotOptions.plotRows = 1000;
	end
	if (plotOptions.plotCols > 1000)
		fprintf('*** Limiting number of plotted columns to 1000\n');
		plotOptions.plotCols = 1000;
	end
	% Pull some options out for easier coding here:
	scalingMainComponent = plotOptions.scalingMainComponent;
	mainSignVectorLength = plotOptions.mainSignVectorLength;
	scalingScdryComponent = plotOptions.scalingScdryComponent;
	gammaPower = plotOptions.gammaPower;

	if (not(ismatrix(localResults)))
		% localResults must be a java matrix
		% Convert the local values back to octave native values, and select
		%  only the rows and cols that will be plotted (gives a big speed up over
		%  coverting all of the values)
		localResultsToPlot = javaMatrixToOctave(localResults, plotOptions.plotStartRow, plotOptions.plotStartCol, ...
				plotOptions.plotRows, plotOptions.plotCols);
		% JIDT jar file is assumed to be on the path since we already have java objects coming in
		mUtils = javaObject('infodynamics.utils.MatrixUtils');
		% Max and min for the entire local profile (use these for the colour box):
		minLocalEntireProfile = mUtils.min(localResults);
		maxLocalEntireProfile = mUtils.max(localResults);
	else
		% Pull out the values we'll be plotting
		localResultsToPlot = localResults(plotOptions.plotStartRow:plotOptions.plotStartRow+plotOptions.plotRows - 1, ...
			plotOptions.plotStartCol:plotOptions.plotStartCol+plotOptions.plotCols - 1);
		% Max and min for the entire local profile (use these for the colour box):
		minLocalEntireProfile = min(min(localResults));
		maxLocalEntireProfile = max(max(localResults));
	end
	% Max and min for what we're plotting (use these for printing out only):
	minLocalOfPlot = min(min(localResultsToPlot));
	maxLocalOfPlot = max(max(localResultsToPlot));
	fprintf('[max,min] for all spacetime local info dynamics profile is [%.3f, %.3f]\n', maxLocalEntireProfile, minLocalEntireProfile);
	fprintf('[max,min] within this particular plot of the profile is    [%.3f, %.3f]\n', maxLocalOfPlot, minLocalOfPlot);
	
	if (plotOptions.scaleColoursToSubsetOfPlot)
		% Scale the colours to the max and min of the points on the plot
		minLocal = minLocalOfPlot;
		maxLocal = maxLocalOfPlot;
	else
		% Scale the colours to the max and min of the entire profile (including points not plotted)
		minLocal = minLocalEntireProfile;
		maxLocal = maxLocalEntireProfile;
	end
	
	if (minLocal < 0)
		% We need a negative component for the colorchart
		if (plotOptions.scaleColoursToExtremes)
			% Put the darkest blues and reds at our max and min respectively
			if (abs(minLocal) > maxLocal)
				% Use a longer length for the red vector than blue
				vNegLength = mainSignVectorLength;
				vPosLength = floor(mainSignVectorLength .* maxLocal ./ abs(minLocal));
			else
				% Use a longer length for the blue vector than red
				vPosLength = mainSignVectorLength;
				vNegLength = floor(mainSignVectorLength .* abs(minLocal) ./ maxLocal);
			end
			bluePosmap = prepareColourmap(vPosLength, true, scalingMainComponent, scalingScdryComponent, gammaPower);
			redNegmap = flipud(prepareColourmap(vNegLength, false, scalingMainComponent, scalingScdryComponent, gammaPower));
			% fprintf('Plotting locals with scaling to extreme values (%d distinct colours for positive, %d for negative)\n', \
			%	vPosLength, vNegLength);
		else
			% Only use the darkest blue/red for which of positive or negative values
			%  had the largest absolute value. For the other, scale the colours
			%  correspondingly.
			% Initialise the full spectrum
			bluePosmap = prepareColourmap(mainSignVectorLength, true, scalingMainComponent, scalingScdryComponent, gammaPower);
			redNegmap = flipud(prepareColourmap(mainSignVectorLength, false, scalingMainComponent, scalingScdryComponent, gammaPower));
			% Then cut down the non-dominant sign's part:
			if (abs(minLocal) > maxLocal)
				% Use a longer length for the red vector than blue; chop blue down
				vPosLength = floor(mainSignVectorLength .* maxLocal ./ abs(minLocal));
				vNegLength = mainSignVectorLength;
				bluePosmap = bluePosmap(1:vPosLength,:);
			else
				% Use a longer length for the blue vector than red; chop red down
				vNegLength = floor(mainSignVectorLength .* abs(minLocal) ./ maxLocal);
				vPosLength = mainSignVectorLength;
				redNegmap = redNegmap(length(redNegmap)-vNegLength + 1:length(redNegmap),:);
			end
			% fprintf('Plotting locals with minor scaled to major (%d distinct colours for positive, %d for negative)\n', \
			%	vPosLength, vNegLength);
		end
		% Construct the colormap with blue for positive and red for negative
		colormap([redNegmap; bluePosmap]);
		% Now, plot the local values with the pre-prepared colormap
		h = imagesc(localResultsToPlot, [minLocal, maxLocal]);
	else
		% We need to pin the minimum value of the blue-only plot to zero.
		bluemap = prepareColourmap(mainSignVectorLength, true, scalingMainComponent, scalingScdryComponent, gammaPower);
		colormap(bluemap);
		% Now, plot the local values with the pre-prepared colormap
		h = imagesc(localResultsToPlot, [0, maxLocal]);
	end
	axis([0.5 (plotOptions.plotCols+0.5) 0.5 (plotOptions.plotRows+0.5)]);
	colorbar
end

