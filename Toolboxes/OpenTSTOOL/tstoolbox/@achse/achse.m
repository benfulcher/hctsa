 function a = achse(argument, varargin)

%tstoolbox/@achse/achse
%   achse class constructor
%
%   Syntax:
%     * a = achse
%       creates default achse object
%     * a = achse(axs)
%       copies achse object axs into a
%     * a = achse(unt)
%       creates achse object using unit unt, with linear spacing, first =
%       0, delta = 1
%     * a = achse(vec)
%       creates achse object with arbitrary spacing, using values in vec
%       as spacing data
%     * a = achse(unt, vec)
%       creates achse object using unit unt with arbitrary spacing, using
%       values in vec as spacing data
%     * a = achse(unt, first, delta)
%       creates achse object with linear spacing, using delta and first
%     * a = achse(unt, first, delta, 'log')
%       creates achse object with logarithmic spacing, using delta and
%       first
%
%   achse used to describe the different dimensions (axes) of a signal
%   object.
%
%   Example:
%     * a = achse(unit('Hz'), 0.01, 10, 'log')
%       creates a logarithmic frequency axis with values 0.01 Hz, 0.1 Hz,
%       1 Hz, 10 Hz
%     * a = achse(label, samplerate)
%       has the same result as
%       a = achse(unit(label), 0, 1/samplerate)
%
%   see also: delta first horzcat label name quantity resolution
%   samplerate scale setname spacing unit
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

a.name = '';					% name of this axis, somtimes equal to its quantity
a.quantity = '';				% the quantity that is denoted by this axis, e.g. Time or Displacement or Profit
a.unit = unit;					% the unit that must be appended to the axis' values
a.resolution = 'linear';		% 'linear' 'logarithmic' 'arbitrary'
a.first = 0;					% determines the starting value
a.delta = 1;					% needed for spacing
a.values = [];					% in case of arbitrary spacing, otherwise empty
a.opt = {}; 					% optional data can be stored here

if nargin == 0	
	a = class(a, 'achse');
elseif isa(argument, 'achse')
	a = argument;
elseif isa(argument, 'double')
	a.resolution = 'arbitrary';
	a.values  = argument(:)';		% ' is needed to ensure that we have a row vector
	a = class(a, 'achse');
elseif isa(argument, 'unit')
	a.unit = argument;
	a.quantity = quantity(argument);
	a.name = a.quantity;
	switch nargin
		case 1
		case 2
			a.resolution = 'arbitrary';
			if isa(varargin{1}, 'double')
				a.values  = varargin{1}(:)'; 	% (:)' creates always a row vector
			else
				warning('need vector of spacing values as second argument')
			end
		case 3
			a.first = varargin{1};
			a.delta = varargin{2};
		case 4
			if strcmp(varargin{3}, 'log') | strcmp(varargin{3}, 'logarithmic')
				a.resolution = 'logarithmic';
				a.first = varargin{1};
				a.delta = varargin{2};
			else
				a.resolution = 'linear';
				a.first = varargin{1};
				a.delta = varargin{2};
			end
		otherwise
			help(mfilename)
			error('wrong number of arguments given')	
	end
	a = class(a, 'achse');
elseif isa(argument, 'char')
	a.unit = unit(argument);
	a.quantity = quantity(a.unit);
	a.name = a.quantity;
	if nargin > 1
		if varargin{1} > 0		% expect double at this point
			a.delta = 1/varargin{1};
		else
			warning('need scalar greater zero as samplerate');
		end
	end
	a = class(a, 'achse');
else
	help(mfilename)
	error('wrong arguments given')
end
