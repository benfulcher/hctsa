function out = EN_entropies(y,whaten,p)
% Computes a range of different types of entropies from the input time series, y
% Uses the wentropy function from Matlab's Wavelet toolbox
% Actually I think this is probably crap
% Ben Fulcher 2009

N = length(y);

switch whaten
	case 'shannon' % Shannon entropy
		out = wentropy(y,'shannon')/N; % scales with N for large N
	case 'logenergy' % Log Energy entropy
		out = wentropy(y,'log energy')/N; % scales with N for large N
    case 'threshold' % Magnitude of the signal greater than some value
        out = wentropy(y,'threshold',p)/N;
    case 'sure'
        out = wentropy(y,'sure',p)/N;
end
    
end