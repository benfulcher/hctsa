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

% function jIntMatrix = octaveToJavaIntMatrix(octaveMatrix)
%
% Convert a native octave matrix to a java int 2D array.
% 
% Assumes the JIDT jar is already on the java classpath - you will get a 
% java classpath error if this is not the case.
%

function jIntMatrix = octaveToJavaIntMatrix(octaveMatrix)

	if (exist ('OCTAVE_VERSION', 'builtin'))
		% We're in octave:
		% Using 'org.octave.Matrix' is much faster than conversion cell by cell
		if ((rows(octaveMatrix)*columns(octaveMatrix)) > 1)
			% Do this the normal way
			tmp = javaObject('infodynamics.utils.OctaveMatrix');
			tmp.loadIntData(reshape(octaveMatrix,1,rows(octaveMatrix)*columns(octaveMatrix)),[rows(octaveMatrix), columns(octaveMatrix)]);
			jIntMatrix = tmp.asIntMatrix();
		else
			% For length 1 arrays, we need to perform a hack here or else
			%  java thinks the length-one array is a scalar.
			% See octaveToJavaDoubleArray for a further description.
			% So instead we'll do this the slow way (doesn't matter for one element only)
			jIntMatrix = javaArray('java.lang.Integer', 1, 1);
			jIntMatrix(1, 1) = octaveMatrix(1);
		end
	else
		% We're in matlab: the native matlab 2D array can be passed to java as is:
		
		jIntMatrix = octaveMatrix;
	end


end

