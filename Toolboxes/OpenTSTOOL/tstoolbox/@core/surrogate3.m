function cout = surrogate3(cin)

%tstoolbox/@core/surrogate3
%   Syntax:
%     * cout = surrogate3(cin)
%
%   Input Arguments:
%     * cin - core object
%
%   create surrogate data for a scalar time series by permuting samples
%   randomly
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


if ndim(cin) > 1
	error('not a scalar time series')
end

y = data(cin);

len = dlens(cin,1);

cout = core(y(randperm(len))); 

