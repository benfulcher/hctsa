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

% function [miHeartToBreath] = runHeartBreathRateKraskovMI()
%
% runHeartBreathRateKraskovMI
% Version 1.0
% Joseph Lizier
% 3/2/2015
%
% Used to explore mutual information in the heart rate / breath rate example of Schreiber --
%  but estimates MI using Kraskov-Stoegbauer-Grassberger estimation.
%
% Note that the paths (to libraries, data, etc) are set assuming this code is run in
% the folder tutorial/sampleExercises/matlabOctave, relative to the main folder of the JIDT distribution.
%
%
% Outputs
% - miHeartToBreath - MI (heart ; breath)


function [miHeartToBreath] = runHeartBreathRateKraskovMI()

	tic;
	
	% Add Octave utilities to the path
	addpath('../../../demos/octave/');

	% Assumes the jar is three levels up - change this if this is not the case
	% Octave is happy to have the path added multiple times; I'm unsure if this is true for matlab
	javaaddpath('../../../infodynamics.jar');

	data = load('../../../demos/data/SFI-heartRate_breathVol_bloodOx.txt');
	
	% Restrict to the samples that Schreiber mentions:
	data = data(2350:3550,:);
	
	% Separate the data from each column:
	heart = data(:,1);
	chestVol = data(:,2);
	bloodOx = data(:,3);
	timeSteps = length(heart);
	
	fprintf('MI for heart rate <-> breath rate for Kraskov estimation with %d samples:\n', timeSteps);

	% Using a KSG estimator for MI is the least biased way to run this:
	miCalc=javaObject('infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov2');

	% Compute an MI value between heart and breath
	
	miCalc.initialise(1,1); % univariate calculation
	miCalc.setProperty('k', '4'); % 4 nearest neighbours for KSG estimator
	miCalc.setObservations(octaveToJavaDoubleArray(heart), ...
				octaveToJavaDoubleArray(chestVol));
	miHeartToBreath = miCalc.computeAverageLocalOfObservations();
	fprintf('MI: = %.3f nats\n', miHeartToBreath);
		
	tElapsed = toc;
	fprintf('Total runtime was %.1f sec\n', tElapsed);
end

