function n = type(d)

%tstoolbox/@description/type
%   return signal type (e.g. 'Correlation function', 'Spectrogram' etc.)
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,1);

n = d.type;
