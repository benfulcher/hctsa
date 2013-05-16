function status = loadall(cwd, filter, handles)

% tstool/loadall
% load all filenames with extensions matching filter in current directory cwd
% the names are inserted into the tstool file list

if ((cwd(end) == '/') | (cwd(end) == '\')) & (length(cwd) > 1)
	cwd = cwd(end-1);
end

lhandle = handles.lboxhandle;

found = 0;
d = dir(cwd);
cwd = fileparts(cwd);
for i=1:length(d)
	filename = d(i).name;
	[path,name,ext,ver] = fileparts(filename);
	if strcmp(cwd, '.')
		fullname = filename;		% keine laestigen Punkte vor Files im lokalen Verzeichnis
	else
		fullname = fullfile(cwd,[name ext]);
	end
	if strcmp(filter, ext) | isempty(filter)
		if exist(fullname, 'file')		
			status = linsert(fullname, handles);
			if status == 0
				found = found + 1;
			end
		end
	end
end

if found == 0
	status = 'no matching files found';
else
	status = 0; 	% OK
end
