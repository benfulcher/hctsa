function [sig]=write(s, filename, varargin)

%tstoolbox/@signal/write
%   Syntax:
%     * write(s, filename) (writes in TSTOOL's own file format)
%     * write(s, filename, 'ASCII')
%     * write(s, filename, 'WAV') (RIFF WAVE FORMAT)
%     * write(s, filename, 'AU') (SUN AUDIO FORMAT)
%     * write(s, filename, 'NLD') (old NLD FORMAT)
%     * write(s, filename, 'SIPP') (si++ file format)
%
%   writes a signal object to file filename (uses matlab's file format)
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt
narginchk(2,3);

if nargin == 2
        [path,nam,ext] = fileparts(filename);
        s=addhistory(s,['write as >' nam '<']);       
	s=setname(s,nam);
 	save(filename, 's');
else
	[path,nam,ext] = fileparts(filename);
	switch varargin{1}
		case {'ASCII', 'ascii', 'ASC', 'asc', 'dat', 'txt'}
		    dat = data(s);
			command = ['save ' filename '  dat -ascii -double'];
			eval(command);
		case {'au', 'AU', 'sun'}
			nd = ndim(s);
			if (nd <= 2)
				dat = data(s);
				dat = dat / (1.01*max(abs(dat(:))));		% avoid clipping
				auwrite(dat, samplerate(s.xaxis{1}), filename);
			else
			 	error('signal is not suitable to write to AU audio format');	   
			end
		case {'wav', 'wave', 'WAVE', 'WAV'}
			nd = ndim(s);
			if (nd <= 2)
				dat = data(s);
				dat = dat / (1.01*max(abs(dat(:))));		% avoid clipping
				wavwrite(dat, samplerate(s.xaxis{1}), filename);
			else
				error('signal is not suitable to write to WAVE audio format');
			end
		case {'MAT', 'mat'}	
			dat = data(s);
			command = ['save ' filename '  dat '];
			eval(command);
		case {'NLD' , 'nld', 'NL' , 'nl'}
			dat = data(s);
			param.yname = yname(s);
			param.yunit = label(yunit(s));
			for i=1:ndim(s)
				a = getaxis(s, i);
				param.xdelta(i) = delta(a);
				param.xfirst(i) = first(a);
				param.xnames{i} = name(a);
				param.xunits{i} = label(a);
			end
			nldwrite(filename, dat, comment(s), param);
		case {'SIPP', 'sipp', 'si', 'SI', 'si++', 'SI++'}
			sippsig = sipp(s);
			write(sippsig, filename);
% 		case {'netCDF', 'nCDF', 'NETCDF', 'netcdf', 'nc' }
% 			error('not yet implemented');
		otherwise	
			error(['cannot write in format ' varargin{1}])
	end
end



sig=s;
