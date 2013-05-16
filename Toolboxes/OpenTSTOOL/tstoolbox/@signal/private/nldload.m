function	[status, data, comment, param] = nldload(filename)

% [status, data, comment, param] = nldload(filename)
% MATLAB-Leseroutine für das NLD-Fileformat
% Christian Merkwirth August 1997
% Matlabs Matrizen (auch mehrdimens.) sind so im Speicher abgelegt, dass
% der Wert von A(3,8,10,12,2345) direkt auf A(2,8,10,12,2345) folgt
% Die Werte des letzten Index liegen am weitesten auseinander

[fid, message] = fopen(filename, 'r', 'ieee-le'); % Open for reading, binary low-endian format


comment = {};
data = [];
param.xdelta = [];			% Schrittweite eines Samplessteps in Einheiten der korresponiderenden xunit, bei 1000 Hz Samplerate = 0.001 s
param.xfirst = []; 		% Lage des ersten Samples auf der Achse, dies muss nicht unbedingt bei Null sein
param.xnames = {};
param.xunits = {};			% n means Samples, other possibilities Hz, s, m, etc.

param.dtype = 'binary';	% 'ascii', 'binary' (double) 'complex' (je 2 Double Values)
param.yname = '';
param.yunit = '';

if fid == -1
	disp(message)
	status = -10;	% load opening error
	return;
else	
	[status, dlens, comment, param, revision] = parseheader(fid);
	if status == 0
		expected = prod(dlens); 	% Erwartete Anzahl Vektoren
		switch param.dtype
			case 'binary'			% erwartet double Werte
				[data,count] = fread(fid, expected, 'float64');	% Einlesen in einen Spaltenvektor
				if ( count < expected )
					warning('not enough samples in body of file, padding with zeros');
					data = [data; zeros(expected - count,1) ];		% if not enough data read, pad with zeros
				end
			case 'complex'
				[data,count] = fread(fid, 2 * expected, 'float64');	% Einlesen in einen Spaltenvektor
				if ( count < 2 * expected )
					warning('not enough samples in body of file, padding with zeros');
					data = [data; zeros(2 * expected - count,1) ];		% if not enough data read, pad with zeros
				end
				data = reshape(data, 2, expected);
				data = data(1,:) + i * data(2,:);
			case 'ascii'
				if length(dlens) <= 2
					[data, count] = fscanf(fid, inf, '%f');
				else
					status = -11;	% ASCII allows max. 2 dimensions
				end
				
		end
		if length(dlens) > 1
			data = reshape(data, dlens);
		end
		fclose(fid);
		status = 0;
	else
		fclose(fid);
		return;
	end				
end	


	

	






