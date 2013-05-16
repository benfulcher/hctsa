function	[status, dlens, comment, param, revision, last_header_block] = parseheader(fid)

% Diese Funktion liest nur die Headerdaten des Files, nicht den Datenteil
% Der Filepointer sollte danach genau auf den Anfang der Daten zeigen
% In Matlab ist die erste Dimension die ist, in der die Daten im Speicher
% am dichtesten liegen. Auf dieser wird per Default auch gearbeitet.
% Diese erste Dimension ist auch die Column-Dimension
% (Columns = Dim. 1, Rows = Dim.2 ...)
% Das gleiche gilt fuer das NLD-Format

comment = {};
param = defaultparam;
dlens = [];		% [3 10 8] means length of first dimension = 3, ...
nd = 0;		% Number of axes (dimensions)
last_header_block  = -1;			% Stichwort bislang noch gar nicht aufgetaucht
line = deblank(fgetl(fid));

if strncmp(line,'P NLD-TSTOOL', 12)		% File-Magic vorhanden ?
	if length(line) > 12
		revision = deblank(line(13:end));	
		revision = fliplr(deblank(fliplr(revision)));	
	else
		revision = '0';
	end
	
	switch revision 
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case '1.0'
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			while ~strcmp(line,'P EOH')
				line = fgetl(fid);
				if feof(fid) ~= 0
					status = -1;		% Fileende erreicht vor dem Magic EOH
					warning('End-of-File reached before Header was finished');
					return;
				end

				rest = deblank(line(2:end));			% nachfolgende Blanks entfernen
				rest = fliplr(deblank(fliplr(rest)));	% fuehrende Blanks entfernen	
				switch line(1)		% erster Buchstabe entscheidet
					case	{'P', 'p'}
						if strcmp(rest,'EOH')		% Header bereits fertig ?
							break
						end
						[token,rem] = strtok(rest, ' :='); % Parametername kann von Blank,= oder : beendet werden
						if ~isempty(rem)
							rem = rem(2:end);		% Blank, = oder : vom Remainder noch entfernen
						end
						
						parsed = 0; 				% Flag wird gesetzt, sobald dieser Parameter behandelt wurde
						if strcmp(token(1),'n')
							length = str2num(rem);
							if length > 0
								nd = nd +1;
								dlens(1, nd) = length;
								param.xdelta(1,nd) = [1];
								param.xfirst(1,nd) = [0];
								param.xnames{nd,1} = token(2:end);
								param.xunits{nd,1} = '';
							end
							parsed = 1;
						else
							token = lower(token);
							if nd > 0	
								if strcmp(token,'xname')
									if isempty(param.xnames{nd})
										param.xnames{nd} = rem;
									else
										warning(['Header Warning : xname already set, skipping new value']);
									end
									parsed = 1;
								end
								if strcmp(token,'xunit')
									if isempty(param.xunits{nd})
										param.xunits{nd} = rem;
									else
										warning(['Header Warning : xunit already set, skipping new value']);
									end
									parsed = 1;
								end
								if strcmp(token,'xfirst')
									param.xfirst(nd) = str2num(rem);
									parsed = 1;
								end
								if strcmp(token,'xdelta')			
									param.xdelta(nd) = str2num(rem);
									parsed = 1;	
								end
							else
								switch token
									case { 'xname', 'xunit', 'xfirst', 'xdelta'}		% diese sollen nicht unter param.extra auftauchen
										parsed = 1;	
								end				
							end
							if strcmp(token,'dtype')
								rem = lower(rem);
								switch rem
									case { 'ascii', 'binary' , 'complex' }
										param.dtype = rem;
									otherwise
										warning(['Data type ' rem ' not supported']);
										status = -2;	% Data type not supported
										return;
								end
								parsed = 1;
							end
							if strcmp(token,'yname')
								if isempty(param.yname)
									param.yname = rem;
								else
									warning(['Header Warning : yname already set, skipping new value']);
								end
								parsed = 1;
							end
							if strcmp(token,'yunit')
								if isempty(param.yunit)
									param.yunit = rem;
								else
									warning(['Header Warning : yunit already set, skipping new value']);
								end
								parsed = 1;
							end
							if strcmp(token,'last_header_block')
								temp = str2num(rem);
								if last_header_block ~= 0 & last_header_block~= 1
									warning('Wrong argument for last_header_block');
								else
									last_header_block = temp;
								end
								parsed = 1;
							end
							if ~parsed
								siz = size(param.extra);
								ende = siz(1);
								param.extra{ende+1,1} = token;
								param.extra{ende+1,2} = rem;
								parsed = 1;
							end
						end
					case 	{'C', 'c'}
						comment{end+1, 1} = rest;
					otherwise
						warning(['Wrong header line : ' line]);
				end 
			end 					% while ~strcmp(line,'P EOH')
			if  nd > 0		
				status = 0;
			else
				warning('No valid axis found');
				status = -3;	% No valid axis found
				return;
			end
			
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
		otherwise
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		    while ~strcmp(line,'P EOH')
            line = fgetl(fid);
            if feof(fid) ~= 0
                    status = -1;            % Fileende erreicht vor dem Magic EOH
                    warning('End-of-File reached before Header was finished');
                    return;
            end
            
            rest = deblank(line(2:end));                    % nachfolgende Blanks entfernen
            rest = fliplr(deblank(fliplr(rest)));   % fuehrende Blanks entfernen    
            switch line(1)          % erster Buchstabe entscheidet
                case    {'P', 'p'}
                        if strcmp(rest,'EOH')           % Header bereits fertig ?
                                break
                        end
                        [token,rem] = strtok(rest, ' :='); % Parametername kann von Blank,= oder : beendet werden
                        if ~isempty(rem)
                                rem = rem(2:end);               % Blank, = oder : vom Remainder noch entfernen
                        end
                        token = lower(token);
                        parsed = 0;                             % Flag wird gesetzt, sobald dieser Parameter behandelt wurde
                        if strcmp(token(1),'n')
                                length = str2num(rem);
                                if length > 0
                                        nd = nd +1;
                                        dlens(1, nd) = length;
                                        param.xdelta(1,nd) = [1];
                                        param.xfirst(1,nd) = [0];
                                        param.xnames{nd,1} = '';
                                        param.xunits{nd,1} = '';
                                end
                                parsed = 1;
                        else
                                if nd > 0       
                                        if strcmp(token,'xname')
                                                if isempty(param.xnames{nd})
                                                        param.xnames{nd} = rem;
                                                else
                                                        warning(['Header Warning : xname already set, skipping new value']);
                                                end
                                                parsed = 1;
                                        end
                                        if strcmp(token,'xunit')
                                                if isempty(param.xunits{nd})
                                                        param.xunits{nd} = rem;
                                                else
                                                        warning(['Header Warning : xunit already set, skipping new value']);
                                                end
                                                parsed = 1;
                                        end
                                        if strcmp(token,'xfirst')
                                                param.xfirst(nd) = str2num(rem);
                                                parsed = 1;
                                        end
                                        if strcmp(token,'xdelta')                       
                                                param.xdelta(nd) = str2num(rem);
                                                parsed = 1;     
                                        end
                                else
                                        switch token
                                                case { 'xname', 'xunit', 'xfirst', 'xdelta'}            % diese sollen nicht unter param.extra auftauchen
                                                        parsed = 1;     
                                        end                             
                                end
                                if strcmp(token,'dtype')
                                        rem = lower(rem);
                                        switch rem
                                                case { 'ascii', 'binary' , 'complex' }
                                                        param.dtype = rem;
                                                otherwise
                                                        warning(['Data type ' rem ' not supported']);
                                                        status = -2;    % Data type not supported
                                                        return;
                                        end
                                        parsed = 1;
                                end
                                if strcmp(token,'yname')
                                        if isempty(param.yname)
                                                param.yname = rem;
                                        else
                                                warning(['Header Warning : yname already set, skipping new value']);
                                        end
                                        parsed = 1;
                                end
                                if strcmp(token,'yunit')
                                        if isempty(param.yunit)
                                                param.yunit = rem;
                                        else
                                                warning(['Header Warning : yunit already set, skipping new value']);
                                        end
                                        parsed = 1;
                                end
                                if ~parsed
                                        siz = size(param.extra);
                                        ende = siz(1);
                                        param.extra{ende+1,1} = token;
                                        param.extra{ende+1,2} = rem;
                                        parsed = 1;
                                end
                        end
                case    {'C', 'c'}
                        comment{end+1, 1} = rest;
                otherwise
                        warning(['Wrong header line : ' line]);
        	end 		% switch line
        end 	% while 
        if  nd > 0  			
        	status = 0;
        else
        	warning('No valid axis found');
        	status = -3;	% No valid axis found
        	return;
        end
	
	end 						% switch revision
else
	warning('Not an NLD File (magic not found)');
	status = -4;	% Not an NLD File (magic not found)
	return;
end
