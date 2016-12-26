function d = setoptparams(d, varargin)

%tstoolbox/@description/setoptparams
%   Syntax:
%     * d = setoptparams(d, nr, param)
%
%   set optional parameter number nr
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,3);


if nargin < 3
	nr = 1;
	param = varargin{1};
else
	nr = varargin{1};
	param = varargin{2};
end

d.optparam{nr} = param;

