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

% function plotRawCa(states, saveIt)
%
% Plot the given raw values of a cellular automata run
%
% Inputs
% - states - 2D array of states of the CA (1st index time goes along the rows, 2nd index cells go across the columns)
% - rule - Rule number (as number or a string) for the rule, used to label the file
% - plotOptions - structure (optional) containing the following variables:
%   - plotRows - how many rows to plot (default is all)
%   - plotCols - how many columns to plot (default is all)
%   - plotStartRow - which row to start plotting from (default is 1)
%   - plotStartCol - which column to start plotting from (default is 1)
% - saveIt - whether to save a file of the image (default false)
% - saveFormat - what file format to save the image in - 'eps' or 'pdf' (default 'eps')

function plotRawCa(states, rule, plotOptions, saveIt, saveFormat)

	if (nargin < 2)
		plotOptions = {};
	end
	if not(isfield(plotOptions, 'plotRows'))
		plotOptions.plotRows = size(states, 1);
	end
	if not(isfield(plotOptions, 'plotCols'))
		plotOptions.plotCols = columns(states);
	end
	if not(isfield(plotOptions, 'plotStartRow'))
		plotOptions.plotStartRow = 1;
	end
	if not(isfield(plotOptions, 'plotStartCol'))
		plotOptions.plotStartCol = 1;
	end
	if (plotOptions.plotRows > 1000)
		fprintf('*** Limiting number of plotted rows to 1000\n');
		plotOptions.plotRows = 1000;
	end
	if (plotOptions.plotCols > 1000)
		fprintf('*** Limiting number of plotted columns to 1000\n');
		plotOptions.plotCols = 1000;
	end
	if (nargin < 3)
		saveIt = false;
	end
	if (nargin < 4)
		saveFormat = 'eps';
	end

	figure(1);
	% set the colormap for black = 1, white = 0
	colormap(1 - gray(2))
	imagesc(states(plotOptions.plotStartRow:plotOptions.plotStartRow+plotOptions.plotRows - 1, ...
		plotOptions.plotStartCol:plotOptions.plotStartCol+plotOptions.plotCols - 1))
	axis([0.5 (plotOptions.plotCols+0.5) 0.5 (plotOptions.plotRows+0.5)]);
	colorbar
	fprintf('Adding colorbar to ensure that the size of the diagram matches that of local info plots\n');
	if (saveIt)
		if (strcmp(saveFormat, 'eps'))
			printDriver = 'epsc'; % to force colour
			fontSize = 32;
		else
			printDriver = saveFormat;
			fontSize = 13;
		end
		set(gca, 'fontsize', fontSize);
		colorbar('fontsize', fontSize);
		if (ischar(rule))
			filename = sprintf('figures/%s-raw.%s', rule, saveFormat);
		else
			filename = sprintf('figures/%d-raw.%s', rule, saveFormat);
		end
		print(filename, sprintf('-d%s', printDriver));
	end
end

