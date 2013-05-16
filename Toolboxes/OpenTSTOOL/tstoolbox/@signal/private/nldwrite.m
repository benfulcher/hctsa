function status = nldwrite(filename, data, comment, param)

% status = nldwrite(filename, data, comment, param)
% data kann auch komplexwertig sein
% MATLAB-Schreibroutine für das NLD-Fileformat
% Christian Merkwirth August 1997

TSTOOL_REVISION = '1.0';		% 30.Sep.1997 cmerk


[fid, message] = fopen(filename, 'w', 'ieee-le'); % Open for reading, binary low-endian format


if fid == -1
	disp(message)
	status = -500;		% Write open error
	return;
else
	fprintf(fid, ['P NLD-TSTOOL ' TSTOOL_REVISION ' \n']);
	
	dims = tsndims(data);
	dlens = tssize(data);
	
	% Parameter des Datenteils schreiben
	if isreal(data)	
		fprintf(fid, ['P dtype binary\n']);
	else
		fprintf(fid, ['P dtype complex\n']);
	end
	fprintf(fid, ['P yname ' param.yname '\n']);	
	fprintf(fid, ['P yunit ' param.yunit '\n']);
	
	% Parameter der einzelnen Dimensionen (Axen) schreiben
	
	for i=1:dims	
		if ~(length(param.xnames) < i)
		    temp1 = param.xnames{i};
		else
		    temp1 = '';	
		end
		temp2 = num2str(dlens(i));
		fprintf(fid, ['P n' temp1 ' ' temp2 ' \n']);
		if ~(length(param.xfirst) < i)
			fprintf(fid, ['P xfirst ' num2str(param.xfirst(i)) '\n']);
		end
		if ~(length(param.xdelta) < i)
			fprintf(fid, ['P xdelta ' num2str(param.xdelta(i)) '\n']);
		end
		if ~(length(param.xunits) < i)
			fprintf(fid, ['P xunit ' param.xunits{i} '\n']);
		end		
	end
		
	% Kommentar anfuegen
	if ~isempty(comment)
		[n,m] = size(comment);	
		for i=1:n
			work = comment{i};		% jeweils eine Zeile
			fprintf(fid, ['C ' work '\n']);	% mit CR
		end
	end
	fprintf(fid, 'P EOH\n');	
	% Header beendet
	
	% Datenteil schreiben
	data = data(:); 	
	if isreal(data)	
		fwrite(fid, data(:), 'double'); 
	else
		rdata = real(data);
		idata = imag(data);
		data = [rdata(:)' ;  idata(:)'];	% komplexe Zahlen stehen wie folgt im File : real imag real imag real imag ...
		fwrite(fid, data(:), 'double'); 
	end								
	status = fclose(fid);
	if status ~= 0
		status = -501	% Write close error
	end
end	


	

	






