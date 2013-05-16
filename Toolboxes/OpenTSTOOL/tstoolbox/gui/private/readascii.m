function l = readascii(filename)

% l = readascii(filename)
% read each line of an ASCII file into list l

l = list;		% create empty list
fid = fopen(filename, 'r');

if fid ~= -1
	while 1
		line = fgetl(fid);
		if ~isstr(line), break, end
			l = append(l, line);
		end
	fclose(fid);
end
