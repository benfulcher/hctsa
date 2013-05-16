function param = defaultparam

% function param = defaultparam
% Erzeugt eine param Struktur


param.xdelta = [];		% Vektor						% Schrittweite eines Samplessteps in Einheiten der korresponiderenden xunit, bei 1000 Hz Samplerate = 0.001 s
param.xfirst = []; 		% Vektor 					% Lage des ersten Samples auf der Achse, dies muss nicht unbedingt bei Null sein
param.xnames = {};		% Cell-Array mit Strings
param.xunits = {};		% Cell-Array mit Strings		% n means Samples, other possibilities Hz, s, m, etc.

param.dtype = 'binary';		% String				% 'ascii', 'binary' (double) 'complex' (je 2 Double Values)
param.yname = '';			% String
param.yunit = '';			% String
param.extra = { };		% Cell-Array mit Strings

