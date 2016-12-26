function [out, T] = BF_pareval(x,y,s,beVocal)
% BF_pareval    Sneaky eval function for parfor loops
%
% parfor loops don't allow the 'eval' function. This sneaky sneaky gets around
% that.
%---INPUTS:
% x and y are possible elements of the string s to be evaluated
% Stores any text output as T

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

if beVocal
    % Any text output from operations is printed to screen
    out = eval(s);
else
    % Text output from operations is suppressed and stored in T
    [T, out] = evalc(s);
end

end
