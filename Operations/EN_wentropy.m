function out = EN_wentropy(y,whaten,p)
% EN_wentropy   Entropy of time series using wavelets.
%
% Uses the wentropy function from Matlab's Wavelet toolbox.
%
%--INPUTS:
% y, the input time series
% whaten, the entropy type:
%               'shannon',
%               'logenergy',
%               'threshold' (with a given threshold),
%               'sure' (with a given parameter).
%               (see the wentropy documentaiton for information)
% p, the additional parameter needed for threshold and sure entropies
%
%---NOTE:
% It seems likely that this implementation of wentropy is nonsense.

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

% ------------------------------------------------------------------------------
%% Check that a Wavelet Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('wavelet_toolbox')

% ------------------------------------------------------------------------------
% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(whaten)
    whaten = 'shannon'; % default
end

N = length(y); % time-series length

switch whaten
	case 'shannon' % Shannon entropy
		out = wentropy(y,'shannon')/N; % scales with N for large N

	case 'logenergy' % Log Energy entropy
		out = wentropy(y,'log energy')/N; % scales with N for large N

    case 'threshold' % Magnitude of the signal greater than some value
        out = wentropy(y,'threshold',p)/N;

    case 'sure'
        % Equivalent to threshold entropy?
        out = wentropy(y,'sure',p)/N;

    otherwise
        error('Unknown entropy type ''%s''', whaten);
end

end
