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
% J.T. Lizier, M. Prokopenko and A.Y. Zomaya, "A framework for the local information dynamics of distributed computation in complex systems", in Guided self-organisation: Inception, ed. M. Prokopenko, pp. 115-158, Springer, Berlin/Heidelberg, 2014; doi: 10.1007/978-3-642-53734-9_5, arXiv:0811.2690.

clear all;

% Set up simulation options:
cells = 10000;
timeSteps = 600;
neighbourhood = 3;
caStates = 2;

% Set up options for information dynamics analysis, and which segment of the CA to plot
measureParams.k=16;
options.saveImages = true;
options.saveImagesFormat = 'pdf';
options.plotOptions.scaleColoursToSubsetOfPlot = true;
scaleColoursToExtremesDefault = false;
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault;
% Turn up the contrast so that the small values aren't disproportionately visible: (0.15, 0.30 was good, except for separable which was better with 0.15, 0.35)
options.plotOptions.scalingMainComponent = 0.15;
options.plotOptions.scalingScdryComponent = 0.30;
options.plotOptions.gammaPower = 0.5;

%%%%%%%%%
% Examining rule 54:
options.plotOptions.plotRows = 35;
options.plotOptions.plotCols = 35;
options.plotOptions.plotStartRow = 150;
options.plotOptions.plotStartCol = 175;
options.seed = 2;
if (exist('initialStates/GsoChapterDemo2013-initialState.txt', 'file'))
	% A file specifying the initial state exists -- this
	%  ensures that Matlab and Octave use the same initial state
	%  (otherwise only Octave recreates the same initial state used in our chapter).
	%  (You can delete/move the initial state file if you want them generated from scratch.)
	options.initialState = load('initialStates/GsoChapterDemo2013-initialState.txt');
elseif (isfield(options, 'initialState'))
	options = rmfield(options, 'initialState');
end
fprintf('\nStarting rule 54 ...\n');
fprintf('\nPlotting active info storage ...\n');
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'active', measureParams, options);
fprintf('\nPress any key when ready for excess entropy ...\n');
pause
options.plotRawCa = false;
measureParams.k=8; % just for excess entropy
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'excess', measureParams, options);
fprintf('\nPress any key when ready for apparent transfer entropy j = 1 ...\n');
pause
% Use the full red scale for transfer and separable info, since we need to see the extreme negative values properly
options.plotOptions.scaleColoursToExtremes = true;
options.plotOptions.scalingScdryComponent = 0.35;
options.plotOptions.scalingMainComponent = 0.35; 
measureParams.k=16; % back to default
measureParams.j = 1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for apparent transfer entropy j = -1 ...\n');
pause
measureParams.j = -1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for the separable information ...\n');
pause
options.plotOptions.scalingMainComponent = 0.15; % Return to previous value
options.plotOptions.scalingScdryComponent = 0.35;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'separable', measureParams, options);
fprintf('\nPress any key when ready to apply to the next rule\n')
pause
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault; % return to default value 
options.plotOptions.scalingScdryComponent = 0.30; % return to previous value

%%%%%%%%%
% Examining rule 110: (60, 60, 220, 230 is good but with only one collision)
options.plotRawCa = true;
options.plotOptions.plotRows = 50;
options.plotOptions.plotCols = 50;
options.plotOptions.plotStartRow = 50+20;
options.plotOptions.plotStartCol = 800+60;
options.seed = 2;
if (exist('initialStates/GsoChapterDemo2013-initialState.txt', 'file'))
	% A file specifying the initial state exists -- this
	%  ensures that Matlab and Octave use the same initial state
	%  (otherwise only Octave recreates the same initial state used in our chapter).
	%  (You can delete/move the initial state file if you want them generated from scratch.)
	options.initialState = load('initialStates/GsoChapterDemo2013-initialState.txt');
elseif (isfield(options, 'initialState'))
	options = rmfield(options, 'initialState');
end
fprintf('\nStarting rule 110 ...\n');
fprintf('\nPlotting active info storage ...\n');
plotLocalInfoMeasureForCA(neighbourhood, caStates, 110, cells, timeSteps, 'active', measureParams, options);
fprintf('\nPress any key when ready for excess entropy ...\n');
pause
options.plotRawCa = false;
measureParams.k=8; % just for excess entropy
plotLocalInfoMeasureForCA(neighbourhood, caStates, 110, cells, timeSteps, 'excess', measureParams, options);
fprintf('\nPress any key when ready for entropy rate ...\n');
pause
measureParams.k=16; % back to default
plotLocalInfoMeasureForCA(neighbourhood, caStates, 110, cells, timeSteps, 'entropyrate', measureParams, options);
fprintf('\nPress any key when ready for apparent transfer entropy j = -1 ...\n');
pause
% Use the full red scale for transfer and separable info, since we need to see the extreme negative values properly
options.plotOptions.scaleColoursToExtremes = true;
measureParams.j = -1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 110, cells, timeSteps, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for the separable information ...\n');
pause
options.plotOptions.scalingScdryComponent = 0.35;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 110, cells, timeSteps, 'separable', measureParams, options);
fprintf('\nPress any key when ready to apply to the next rule\n')
pause
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault; % return to default value 
options.plotOptions.scalingScdryComponent = 0.30; % return to previous value

%%%%%%%%%
% Examining rule phi_par:
% (75, 75, 12, 300 - is ok, but another collision would be good)
% (100, 100, 12, 500 - is very interesting, but perhaps slightly too complicated. Dynamics will lead to some strange colouring also)
options.plotRawCa = true;
options.plotOptions.plotRows = 70;
options.plotOptions.plotCols = 70;
options.plotOptions.plotStartRow = 17;
options.plotOptions.plotStartCol = 585;
options.seed = 2;
if (exist('initialStates/GsoChapterDemo2013-initialState-phipar.txt', 'file'))
	% A file specifying the initial state exists -- this
	%  ensures that Matlab and Octave use the same initial state
	%  (otherwise only Octave recreates the same initial state used in our chapter).
	%  (You can delete/move the initial state file if you want them generated from scratch.)
	options.initialState = load('initialStates/GsoChapterDemo2013-initialState-phipar.txt');
elseif (isfield(options, 'initialState'))
	options = rmfield(options, 'initialState');
end
phi_par = 'feedffdec1aaeec0eef000a0e1a020a0';
phi_par_neighbourhood = 7;
phi_par_cells = 30000;
phi_par_timeSteps = 200;
measureParams.k=10; % Shorter for phi_par
fprintf('\nStarting rule phi_par ...\n');
fprintf('\nPlotting active info storage ...\n');
plotLocalInfoMeasureForCA(phi_par_neighbourhood, caStates, phi_par, phi_par_cells, phi_par_timeSteps, 'active', measureParams, options);
options.plotRawCa = false;
fprintf('\nPress any key when ready for apparent transfer entropy j = -1 ...\n');
pause
% Use the full red scale for transfer and separable info, since we need to see the extreme negative values properly
options.plotOptions.scaleColoursToExtremes = true;
measureParams.j = -1;
plotLocalInfoMeasureForCA(phi_par_neighbourhood, caStates, phi_par, phi_par_cells, phi_par_timeSteps, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for complete transfer entropy j = -1 ...\n');
pause
plotLocalInfoMeasureForCA(phi_par_neighbourhood, caStates, phi_par, phi_par_cells, phi_par_timeSteps, 'transfercomplete', measureParams, options);
fprintf('\nPress any key when ready for apparent transfer entropy j = -3 ...\n');
pause
measureParams.j = -3;
plotLocalInfoMeasureForCA(phi_par_neighbourhood, caStates, phi_par, phi_par_cells, phi_par_timeSteps, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for the separable information ...\n');
pause
options.plotOptions.scalingMainComponent = 0.35; % Make it easier to see stronger negatives
options.plotOptions.scalingScdryComponent = 0.45;
plotLocalInfoMeasureForCA(phi_par_neighbourhood, caStates, phi_par, phi_par_cells, phi_par_timeSteps, 'separable', measureParams, options);
fprintf('\nPress any key when ready to apply to the next rule\n')
pause
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault; % return to default value 
options.plotOptions.scalingScdryComponent = 0.30; % return to previous value
options.plotOptions.scalingMainComponent = 0.15; % Return to previous value
measureParams.k=16; % back to default

%%%%%%%%%
% Examining rule 22:
options.plotOptions.plotRows = 50;
options.plotOptions.plotCols = 50;
options.plotOptions.plotStartRow = 150;
options.plotOptions.plotStartCol = 175;
options.seed = 2;
if (exist('initialStates/GsoChapterDemo2013-initialState.txt', 'file'))
	% A file specifying the initial state exists -- this
	%  ensures that Matlab and Octave use the same initial state
	%  (otherwise only Octave recreates the same initial state used in our chapter).
	%  (You can delete/move the initial state file if you want them generated from scratch.)
	options.initialState = load('initialStates/GsoChapterDemo2013-initialState.txt');
elseif (isfield(options, 'initialState'))
	options = rmfield(options, 'initialState');
end
fprintf('\nStarting rule 22 ...\n');
fprintf('\nPlotting active info storage ...\n');
options.plotRawCa = true;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 22, cells, timeSteps, 'active', measureParams, options);
fprintf('\nPress any key when ready for excess entropy ...\n');
pause
options.plotRawCa = false;
measureParams.k=8; % just for excess entropy
plotLocalInfoMeasureForCA(neighbourhood, caStates, 22, cells, timeSteps, 'excess', measureParams, options);
fprintf('\nPress any key when ready for entropy rate ...\n');
pause
measureParams.k=16; % back to default
plotLocalInfoMeasureForCA(neighbourhood, caStates, 22, cells, timeSteps, 'entropyrate', measureParams, options);
fprintf('\nPress any key when ready for apparent transfer entropy j = 1 ...\n');
pause
% Use the full red scale for transfer and separable info, since we need to see the extreme negative values properly
options.plotOptions.scaleColoursToExtremes = true;
measureParams.j = 1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 22, cells, timeSteps, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for the separable information ...\n');
pause
options.plotOptions.scalingScdryComponent = 0.35;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 22, cells, timeSteps, 'separable', measureParams, options);
fprintf('\nPress any key when ready to apply to the next rule\n')
pause
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault; % return to default value 
options.plotOptions.scalingScdryComponent = 0.30; % return to previous value

