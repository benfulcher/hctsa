function rs = level_adaption(s, timeconstants, dynamic_limit, threshold)

%tstoolbox/@signal/level_adaption
%   Syntax:
%     * level_adaption(s, timeconstants, dynamic_limit, threshold)
%
%   Each channel of signal s is independently divided by a scaling factor
%   that adapts to the current level of the samples in this channel. The
%   adaption process is simulated using a cascade of feedback loops
%   (P�schel 1998) which consists of low pass filters with time constants
%   given as second argument to this function. The number of time
%   constants given determines the number of feedback loops that are used.
%
%   Higher values for time constants will result in slower adaption speed.
%   Short time changes in the signal will be transmitted almost linearily.
%   In each feedback loop, a nonlinear compressing characteristic (see
%   Stefan M�nkner 1993) limits the signal values to be within
%   [-dynamic_limit dynamic_limit]. A low value for dynamic_limit will
%   introduce nonlinear distortions to the signal.
%
%   To prevent the feedback loops from adapting to a zero level (in case
%   all input values are zero), a tiny threshold is given as 4th argument.
%   The scaling factors will not shrink below this threshold.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,4)

if nargin < 4
	threshold = 0.0000001;
end
if nargin < 3
	dynamic_limit = 10;
end

rs = signal(core(level_adaption(data(s), timeconstants, dynamic_limit, threshold)), s);	
rs = addhistory(rs,  ['Level adaption']);
rs = addcommandlines(rs, 's = level_adaption(s', timeconstants, dynamic_limit, threshold);

