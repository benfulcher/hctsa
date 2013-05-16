function  status = writeascii(filename, l)

% status = writeascii(filename, l)
% write each line of list l into an ASCII file filename

status = 0;
fid = fopen(filename, 'w');

if fid ~= -1
	for i=1:length(l)
		line = get(l,i);
		fprintf(fid, ['%s \n'], line);
	end
	fclose(fid);
else
	status = -1;
end
