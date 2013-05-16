function cout = diff(cin, nth, delta)

%tstoolbox/@core/diff
%   Syntax:
%     * cout = diff(cin, nth, delta)
%
%   Input Arguments:
%     * cin - core object
%     * nth - number of derivations
%     * delta - time difference between to signal values
%
%   nth numerical derivative along dimension 1 when data was sampled
%   equidistantly with samplerate = 1/delta
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

cout = core(diff(data(cin),nth,1)/(delta^nth));

