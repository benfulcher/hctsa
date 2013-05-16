function cout = int(cin, delta)

%tstoolbox/@core/int
%   Syntax:
%     * cout = int(cin, delta)
%
%   Input Arguments:
%     * cin - core object
%     * delta - time period between two data samples
%
%   numerical integration along dimension 1 when data was sampled
%   equidistantly with samplerate = 1/delta
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

tmp = cumsum(data(cin),1)*delta;
tmp = [zeros(1, dlens(cin,2)) ; tmp]; 
cout = core(tmp);
