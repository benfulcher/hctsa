function d = description(varargin)

%tstoolbox/@description/description
%   description class constructor Syntax:
%     * d = description()
%     * d = description(name)
%     * d = description(name, type)
%     * d = description(yunit)
%
%   An object of class description contains auxiliary descriptive
%   information for a signal, e.g. information about data unit, creator,
%   how the signal should be plotted, a user specified comment text, a
%   processing history and the commandlines that were used to generate
%   this signal
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


d.label = '';				% may be used to freely by the user, will be continued when producing a new signal (optional)
d.name = '';			  	% set when loaded from file, otherwise empty, will not be continued in processing (optional)
d.type = '';		  		% 'Spectrum', 'Time-Signal', 'Matrix', 'Image' , 'Spektrogramm' etc. (optional)
d.plothint = '';			% Recommended way of plotting when displayed with view (optional)
d.comment = list;
d.history = list;
d.creator = '';
d.yname = '';
d.yunit = unit;
d.commandlines = list;
d.optparam = {};		  	% optional data can be stored here

d = class(d, 'description');

switch nargin
	case 1 
		if isa(varargin{1}, 'char')
			d = setname(d, varargin{1});
		elseif isa(varargin{1}, 'unit')
			d = setyunit(d, varargin{1});
		else 
			error('Wrong type of argument(s) given');
		end
	case 2
		d = setname(d, varargin{1});
		d = settype(d, varargin{2});
end















