function nname = newname(oldname, praefix, suffix, newpath, newext)

% returns a filename that does not yet exist using the old name and extension, 
% a praefix and a suffix and if necessary a number

[path,name,ext,ver] = fileparts(oldname);

if nargin < 1, help(mfilename); return; end

if nargin == 1
	praefix = '';
	suffix  = '';
elseif nargin == 2
	suffix = '';
elseif nargin == 4
	path = newpath;
elseif nargin == 5
	ext = newext;
	path = newpath;
end

nname = fullfile(path, [praefix name suffix ext]);

i = 0;
while exist(nname, 'file')
	i = i+1;
	nname = fullfile(path, [praefix name suffix '_' int2str(i) ext]);
end
