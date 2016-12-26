function param = optparams(d, varargin)

%tstoolbox/@description/optparams
%   Syntax:
%     * param = optparams(d, nr)
%
%   get optional parameter number nr
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,2);

if isempty(d.optparam)
	param = [];
	return
end

if nargin < 2
	nr = 1;
else
	nr = varargin{1};
end

param = d.optparam{nr}; 

