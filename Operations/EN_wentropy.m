% EN_wentropy
% 
% Uses the wentropy function from Matlab's Wavelet toolbox to
% estimate the entropy of the input time series.
% 
% INPUTS:
% y, the input time series
% whaten, the entropy type:
%               'shannon',
%               'logenergy',
%               'threshold' (with a given threshold),
%               'sure' (with a given parameter).
%               (see the wentropy documentaiton for information)
% p, the additional parameter needed for threshold and sure entropies
% 
% Author's cautionary note: it seems likely that this implementation of wentropy is nonsense.

function out = EN_wentropy(y,whaten,p)
% Ben Fulcher, 2009

% Check inputs
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
        out = wentropy(y,'sure',p)/N;
    otherwise
        error('Unknown entropy type ''%s''', whaten);
end
    
end