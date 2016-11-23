function d = data(c, varargin)

%tstoolbox/@core/data
%   Syntax:
%     * d = data(c, varargin)
%     * c=core object
%
%   Input Arguments:
%     * varargin - selector string for data-elements in matlab notation
%
%   Return signal's data values
%   With no extra arguments, data returns the data array of a signal
%   object
%   Another possible call is : data(signal, ':,:,1:20')
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,2);

if nargin == 1
	d = c.data;
elseif isa(varargin{1}, 'char') 
	eval(['d = c.data(' varargin{1} ')']);
else
	error('Wrong type of argument(s) given');
end
