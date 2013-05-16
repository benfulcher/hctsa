function s = signal(argument, varargin)

%tstoolbox/@signal/signal
%   Syntax:
%     * s = signal(array)
%       creates a new signal object from a data array array the data
%       inside the object can be retrieved with x = data(s);
%     * s = signal(array, achse1, achse2, ...)
%       creates a new signal object from a data array array, using achse1
%       etc. as xachse entries
%     * s = signal(array, unit1, unit2, ...)
%       creates a new signal object from a data array 'array', using unit1
%       etc. to create xachse objects
%     * s = signal(array, samplerate1, samplerate2, ...)
%       creates a new signal object from a data array array, using as
%       xunit 's' (second) and scalar samplerate1 as samplerate(s)
%
%   A signal object contains signal data, that is a collection of real or
%   complex valued samples. A signal can be one or multi-dimensional. The
%   number of dimensions is the number of axes that are needed to describe
%   the the data.
%
%   An example for an one-dimensional signal is a one-channel measurement
%   (timeseries), or the power spectrum of a one-channel measurement. An
%   example for a two-dimensional signal is a twelve-channel measurement,
%   with one time axis and a 'channel' axis. Another example for a
%   two-dimensional signal is a short time spectrogramm of a time series,
%   where we have a time axis and a frequency axis.
%
%   Each axis can have a physical unit(e.g. 's' or 'Hz'), a starting point
%   and a step value. E.g. if a time-series is sampled with 1000 Hz,
%   beginning at 1 min 12 sec, the unit is 's', the starting point is 72
%   and the step value (delta) is 0.001.
%
%   But not only the axes have physical units, also the sample value
%   themselve can have a unit, maybe 'V' or 'Pa', depending on what the
%   sampled data represent (=> yunit)
%
%   All units are stored as objects of class 'unit', all axes are stored
%   as objects of class 'achse' (this somewhat peculiar name was chosen
%   because of conflicts with reserved matlab keywords 'axis' and 'axes',
%   which otherwise would have been the first choice).
%
%   Example for creating a 2-dimensional signal with y-unit set to 'Volt',
%   the first dimension's unit is 'second' (time), the second dimension's
%   unit is 'n' (Channels).
%
%   Examples:
%     *
%tmp = rand(100, 10);
%
%s = signal(tmp, unit('s'), unit('n'));
%s = setyunit(s, unit('V'));
%s = addcomment(s, 'Example signal with two dimensions')
%
%     * Loading from disk
% s = signal(filename)
%       loads a previously stored signal object
%     * Importing from other file formats:
%ASCII: s = signal('data/spalte1.dat', 'ASCII')
%WAVE: s = signal('data/Sounds/hat.wav', 'WAVE')
%AU (SUN AUDIO): s = signal('data/Sounds/hat.au', 'AU')
%(old) NLD-Format : s = signal('test.nld', 'NLD')
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

if nargin < 1
  argument=zeros(16,1);                 % generating dummy signal,
                                        % not very useful, but
                                        % imported
end
  
  
if isa(argument, 'char')		% loading file with name given by argument
	if nargin == 1
		load(argument, '-mat');
		[path,nam,ext,ver] = fileparts(argument);
		% if isempty(name(s))			% give only new name if old one is empty
		s = setname(s, nam);			
		%end
	else
		[path,nam,ext,ver] = fileparts(argument);
		switch varargin{1}
			case {'ASCII', 'ascii', 'ASC', 'asc', 'dat', 'txt'}
				load(argument,'-ascii'); 
				eval(['dat = ' nam ';']);
				[s,c,d] = create(dat);	s = class(s, 'signal', c,d);
				s = addhistory(s,  ['Imported from ASCII file ''' argument '''']);
			case {'au', 'AU', 'sun'}
				[siz, Fs, bits] = auread(argument,'size');
				samples = siz(1);
				channels = siz(2);
				dat = auread(argument);
				[s,c,d] = create(dat);	s = class(s, 'signal', c,d);
				if Fs > 0
					s = setaxis(s, 1, achse(unit('s'), 0, 1/Fs));		% give signal time axis
				end
				s = addhistory(s, { ['Imported from SUN audio file ''' argument ''''] ...
			    					['Number of channels : ' num2str(channels)] ...
			    					['Bits per Sample : ' num2str(bits)] });
			case {'wav', 'wave', 'WAVE', 'WAV'} 	% Laden von Sound-Dateien im WAVE-Format
				[siz, Fs, bits] = wavread(argument,'size');
				samples = siz(1);
				channels = siz(2);
				dat = wavread(argument);
				[s,c,d] = create(dat);	s = class(s, 'signal', c,d);
				if Fs > 0
					s = setaxis(s, 1, achse(unit('s'), 0, 1/Fs));		% give signal time axis
				end
				s = addhistory(s, { ['Imported from WAVE file ''' argument ''''] ...
									['Number of channels : ' num2str(channels)] ...
									['Bits per Sample : ' num2str(bits)] });
			case {'NLD', 'nld' }		% Laden des alten NLD-Formates
				[status, dat, comment, param] = nldload(argument);
				if status==0 % SUCCESS
					[s,c,d] = create(dat);	s = class(s, 'signal', c,d);
					for i=1:ndim(s)
						a = achse(unit(param.xunits{i}), param.xfirst(i), param.xdelta(i));
						a = setname(a, param.xnames{i}); 
						s = setaxis(s, 1, a);
					end
					s = setyname(s, param.yname);
					s = setyunit(s, param.yunit);
					s = addcomment(s, comment);
					s = addhistory(s, ['Imported from NLD file ''' argument '''']);
				else
					error(['could not read file ' argument ' as NLD-file']);
				end	
                        case {'netCDF', 'nCDF', 'NETCDF', 'netcdf', 'nc' }
% Check for installed netCDF tools from
% http://www.marine.csiro.au/sw/matlab-netcdf.html
                                 if (exist('mexcdf53','file')==3) & ... 
 				      (exist('getnc','file')==2),
 				  dat = getnc(argument); 
 				  [s,c,d] = create(dat);	s = class(s, 'signal', c,d);
 				  s = addhistory(s,  ['Imported from netCDF file ''' argument '''']);
 				else
 				  error('no valid netCDF interface installed');
 				end;
			otherwise
				error('wrong type of arguments')
		end
%%%% <--- netCDF-Patch-new ---end


		s = setname(s, nam);
	end
elseif isa(argument, 'double')		% Erzeugen aus einem einfachen multidimensionalen Array
	hstring = 'Imported from MATLAB workspace'; 		% text to append to history
	[s,c,d] = create(argument); s = class(s, 'signal', c,d);
	nd = ndim(s);
	for i=1:nd
		if (nargin > i) & (isa(varargin{i}, 'achse'))   % construct signal's achse form achse objects
			if strcmp(resolution(varargin{i}), 'arbitrary')
				if length(spacing(varargin{i})) < dlens(s,i)
					error('not enough axis spacing values for this data array');
				else
					s = setaxis(s, i,varargin{i});	
				end
			else
				s = setaxis(s, i,varargin{i});	
			end
		elseif (nargin > i) & (isa(varargin{i}, 'unit'))
			s = setaxis(s, i,achse(varargin{i}));		
		elseif (nargin > i) & (isa(varargin{i}, 'double'))  % construct signal's axis with specified samplerate
			s = setaxis(s, i,achse(unit('s'), 0, 1/varargin{i}));		
		elseif (nargin > i) & (isa(varargin{i}, 'char'))	% use custom specified string for history
			hstring = varargin{i};
		end
	end
	s = addhistory(s, hstring);
elseif isa(argument, 'core')		
	% s = signal( <new core object>, <old signal object>)
	% Dieser Aufruf ist den Arbeitsroutinen dieser Klasse vorbehalten !!!!
	% und nicht fuer den Benutzer von der Kommandozeile aus gedacht, 
	% daher ist er auch nicht in der help dokumnetiert.
	% Es kann u.U. ein signal Objekt erzeugt werden, das nicht vollstaendig initialisiert ist ! 
	% Dies bleibt der aufrufenden Routine ueberlassen
	% expect core and signal object as input arguments
	
	s.xaxes = varargin{1}.xaxes(1:min(ndim(varargin{1}), ndim(argument)));	% do not copy all xaxes
	d = varargin{1}.description;
	d = setname(d, ''); 	% don't continue name 	
	if strcmp(class(optparams(d, 1)), 'struct')
		d = setoptparams(d, 1, []); % don't continue optional parameters (cmerk 3. Aug. 1999)
	end
	
	d = setoptparams(d, {}); % don't continue optparam	
	d = setplothint(d, ''); % don't continue plothint 	
	s = class(s, 'signal', argument , d);
else
	error('wrong type of arguments')
end

function [s,c,d] = create(dat)
s.xaxes = {};					% cell array of achses for each dimension of the signal
d = description;				% produce an empty description object
c = core(dat);					% produce an signal core with the data given by argument
nd = ndim(c);
for i=1:nd						% create all xaxes with default achse object
	s = addaxis(s);
end
