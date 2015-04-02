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

% To recreate the plots in:
% J.T. Lizier and J.R. Mahoney, "Moving frames of reference, relativity and invariance in transfer entropy and information dynamics", Entropy, vol. 15, no. 1, p. 177-197, 2013; doi: 10.3390/e15010177

clear all;

% Set up options for TE j=-1 cell to the right, and which segment of the CA to plot
measureParams.k=16;
measureParams.j = -1;
options.plotOptions.plotRows = 60;
options.plotOptions.plotCols = 60;
options.plotOptions.plotStartRow = 125;
options.plotOptions.plotStartCol = 125;
options.seed = 1;
if (exist('initialStates/MovingFrameDemo2013-initialState.txt', 'file'))
	% A file specifying the initial state exists -- this
	%  ensures that Matlab and Octave use the same initial state
	%  (otherwise only Octave recreates the same initial state used in our chapter).
	%  (You can delete/move the initial state file if you want them generated from scratch.)
	options.initialState = load('initialStates/MovingFrameDemo2013-initialState.txt');
elseif (isfield(options, 'initialState'))
	options = rmfield(options, 'initialState');
end
options.saveImages = true;
options.plotOptions.scaleColoursToSubsetOfPlot = true;
% Turn up the contrast so that the small values aren't disproportionately visible:
options.plotOptions.scalingMainComponent = 0.1;
options.plotOptions.scalingScdryComponent = 0.15;
options.movingFrameSpeed = 0;

% Examining rule 54, stationary frame:
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'active', measureParams, options);
fprintf('\nPress any key when ready for apparent transfer entropy ...\n');
pause
measureParams.j = -1;
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for complete transfer entropy ...\n');
pause
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfercomplete', measureParams, options);
fprintf('\nPress any key when ready for apparent transfer entropy j=1 channel...\n');
pause
measureParams.j = 1;
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for complete transfer entropy ...\n');
pause
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfercomplete', measureParams, options);
fprintf('\nCopy the figures from the figures directory, then press any key')
fprintf('\nwhen ready for the measures in a moving frame of reference\n')
pause

options.movingFrameSpeed = 1;
% Examining rule 54, f=1 frame:
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'active', measureParams, options);
fprintf('\nPress any key when ready for apparent transfer entropy ...\n');
pause
measureParams.j = -1;
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for complete transfer entropy ...\n');
pause
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfercomplete', measureParams, options);
% Also plot TE for channel j=0 in the moving frame:
measureParams.j = 0;
fprintf('\nPress any key when ready for apparent transfer entropy with j=0...\n');
pause
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for complete transfer entropy with j=0...\n');
pause
plotLocalInfoMeasureForCA(3, 2, 54, 10000, 600, 'transfercomplete', measureParams, options);

