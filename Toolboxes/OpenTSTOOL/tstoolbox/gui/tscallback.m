function tscallback(action, varargin)


% tstool/tscallback
% callback handler for tstool graphical user interface
% action is a string that specifies the action to be done
% should only be invoked by the tstool figure

if nargin < 1, help(mfilename), return, end 

try
	global TSTOOLdatapath 

	
	fighandle = findobj('Tag', 'TSTOOL');  % Filename ist in den UserData der TSTOOL-Figure gespeichert
	handles = get(fighandle, 'UserData');	% structure with handles to the different menus/controls of the figure

	fighandle 		= handles.fighandle;				
	lboxhandle		= handles.lboxhandle;
	loadhandle		= handles.loadhandle;
	lallhandle		= handles.lallhandle;
	recopthandle 		= handles.recopthandle;
	plotmodehandle = handles.plotmodehandle;
	defwindowhandle = handles.defwindowhandle;

	answer=get(loadhandle,'UserData');
	TSTOOLfilter=answer{1};
	cwd=answer{2};
	datafiles=get(lboxhandle,'UserData');
%	currfilehandle = handles.currfilehandle;

	% Mehrere UserData Fields des TSTOOLs sind belegt :
	% der TSTOOL-Figure selbst, enthaelt eine Struktur mit Handle-Nummern der verschiedenen Menus/Controls des tstool
	% der Listbox => momentan selektierter Filename
	% des Load-Menus => das momentane Suchverzeichnis (.) und Schreibverzechnis (.)
	% des Load-All Menus => die momentane Filterendung für Filenamen (.sgo)
	% des Embedding-Options-Menus => Einbettungsparameter (Vektor aus dim und delay)

	% FIXME cmerk 27.11.1998 bitte aendern auf neue Buttontypes ('Checked' == 'on/'off')
	% des Plot-Menus => Flag, ob das Ergebnis jeder Aktion automatisch geplottet werden soll
	% der Plot-Axis => preview Modus, wenn dort nicht 'small' steht, wird Plot-Ausgabe in eigenes Fenster geleitet

	% Setzen einiger Defaultwerte

%	filter = '.sgo';			% Extension, an der ein TSTOOL-File erkannt wird	
	dim = 3; delay = 1; shift = 1; recwindowtype = 'Rect';		% Einbettungsparameter
	global ploteverytime plotmode;
	fileselboxXY = [300 300];		% Defaultkoordinaten der Fileselectbox
	ScreenSize = get(0, 'ScreenSize');
	fileselboxXY = [ceil(ScreenSize(3)/4), ceil(ScreenSize(4)/4)];

	%Settings = get(0, 'UserData');
	%TSTOOLpath = Settings{1};
	[TSTOOLpath,dummy,dummy,dummy]=fileparts(which('units.mat'));

%	global cwd % Aktuelles Suchverzeichnis und aktuelles Schreibverzechnis
%	filter = get(lallhandle, 'UserData');

	% Parameters for time-delay reconstruction
%	embparams = get(recopthandle, 'UserData'); 
%	dim = str2num(embparams{1}); delay = str2num(embparams{2});
%	shift = str2num(embparams{3}); recwindowtype = embparams{4};

	% Default window type for spectral applications (not identical to windowtype for reconstruction)
	windowtype = get(defwindowhandle, 'UserData');

	ploteverytime = 0;
	if strcmp(get(plotmodehandle(1), 'Checked'), 'on')
		plotmode = 'small';
		ploteverytime = 1;
	end
	
	if strcmp(get(plotmodehandle(2), 'Checked'), 'on')
		plotmode = 'large';
		ploteverytime = 1;
	end

	% Aktionen werden unterteilt in solche, die einen unter UserData abgelegten Filenamen
	% benötigen (wie die meisten Rechenmodule, und solche, die den Filenamen nicht benötigen,
	% weil sie z.B. diesen erst setzen (wie load oder select)

	% In folgendem Abschnitt werden Callbacks behandelt, die keinen Filenamen benötigen

	switch action
	  case 'exit'
	    savesettings(handles);
	    close(gcf);
	  case 'load'
	  [fname, pname] = uigetfile(fullfile(cwd,['*' TSTOOLfilter]), ...
				     'Load data file', fileselboxXY(1), ...
				     fileselboxXY(2));

	  if fname ~= 0
	    fullname = fullfile(pname, fname);
	    if exist(fullname, 'file')	
	      cwd=pname;
	      set(loadhandle, 'UserData', {TSTOOLfilter cwd});
	      s=signal(fullname);
	      
	      [status,datafiles]=siglevel(s,datafiles);
	      datafiles=sortdatafiles(datafiles);
	      filllistbox(datafiles,handles,name(s));
	      write(s,newname(fname,'','',TSTOOLdatapath,TSTOOLfilter));
	      % datafiles=lload(fullname, handles,datafiles);
	      if ploteverytime, tsplot([pname fname], plotmode); end
	    else
	      warndlg(['Error finding file ' fullname]);
	    end
	  end
	  
	 case 'loadall'
	  % title = 'Load mutiple files';
	  % lineNo = 1;
	  % prompt = {'Filter'};
	  % defaults = {fullfile(cwd, ['*' filter])};
	  % answer = inputdlg(prompt,title,lineNo,defaults);
	  % if ~isempty(answer)
	  [status,datafiles] = ...
	      loadallfiles(fullfile(TSTOOLdatapath, ...
				    ['*' TSTOOLfilter]), TSTOOLfilter);
	  cla
	  if status ~= 0
	    warndlg(status);
	    set(handles.lboxhandle,'String',{});
	    set(handles.lboxhandle,'Value',0);
	    set(handles.lboxhandle,'UserData',{});
	  else
	    datafiles=sortdatafiles(datafiles);
	    filllistbox(datafiles,handles,1);
	  end
	  % end
	 case 'forall'       % apply macro to all files in listbox
	  if exist(fullfile(TSTOOLdatapath, 'scripts', 'macro.m'))
	    clear('macro'); 	% make shure newest macro.m is executed
	    String = get(lboxhandle, 'String');
	    
	    if  ~isempty(String)  
	      for i=1:length(String)
            	currentfile = getcurrentfile(TSTOOLdatapath,String{i});
            	eval('through(''macro'', currentfile);', ''); % try
                                                              % to
                                                              % catch
                                                              % errors
	      end   
	    else
	      warndlg('No signals in list');
	    end    		
	  else
	    warndlg('No macro defined yet');
	  end
	 case 'select'		% Eine anderer Eintrag der
				% Listbox wurde angewaehlt
	   Value = get(lboxhandle, 'Value');
	   String = get(lboxhandle, 'String');
	   if ~isempty(String)  
	     if ploteverytime
	       tsplot(fullfile(TSTOOLdatapath, [getlastentry(datafiles,Value) TSTOOLfilter]), 'small');
	     end
	       %setcurrentfile(String{Value}, handles.lboxhandle, handles.currfilehandle);
	       %set(lboxhandle, 'UserData', String{Value} );
	     else
	       %setcurrentfile('', handles.lboxhandle, handles.currfilehandle);
	       %set(lboxhandle, 'UserData', '');
	       cla
	   end
	 
	 
	 case 'impascii'
	  [fname, pname] = uigetfile( fullfile(cwd, '*.dat'), 'Load ASCII file', fileselboxXY(1), fileselboxXY(2));
	  if fname ~= 0
	    cwd=pname;
	    set(loadhandle, 'UserData', {TSTOOLfilter cwd});
	    sig = signal([pname fname], 'asc');
	    newfilename = newname([pname fname], '' , '', TSTOOLdatapath, TSTOOLfilter);
	    sig=write(sig, newfilename);
	    if ~isempty(newfilename)
	      datafiles=linsert(newfilename, handles,datafiles,sig,0);
	      if ploteverytime, tsplot(newfilename, plotmode); end
	    end
	  end
	 case 'impnld'
	  [fname, pname] = uigetfile( fullfile(cwd, '*.nld'), 'Load NLD file', fileselboxXY(1), fileselboxXY(2));
	  if fname ~= 0
	    cwd=pname;
	    set(loadhandle, 'UserData', {TSTOOLfilter cwd});
	    sig = signal([pname fname], 'NLD');
	    newfilename = newname([pname fname], '' , '', TSTOOLdatapath, TSTOOLfilter);
	    sig=write(sig, newfilename);
	    if ~isempty(newfilename)
	      datafiles=linsert(newfilename, handles,datafiles,sig,0);
	      if ploteverytime, tsplot(newfilename, plotmode); end
	    end
	  end
	 case 'impwav'
	  [fname, pname] = uigetfile( fullfile(cwd, '*.wav'), 'Load wave file', fileselboxXY(1), fileselboxXY(2));
	  if fname ~= 0
	    cwd=pname;
	    set(loadhandle, 'UserData', {TSTOOLfilter cwd});
	    sig = signal([pname fname], 'wav');
	    newfilename = newname([pname fname], '' , '', TSTOOLdatapath, TSTOOLfilter);
	    sig=write(sig, newfilename);
	    if ~isempty(newfilename)
	      datafiles=linsert(newfilename, handles,datafiles,sig,0);
	      if ploteverytime, tsplot(newfilename, plotmode); end
	    end
	  end
	 case 'impau'
	  [fname, pname] = uigetfile( fullfile(cwd, '*.au'), 'Load SUN audio file', fileselboxXY(1), fileselboxXY(2));
	  if fname ~= 0
	    sig = signal([pname fname], 'au');
	    newfilename = newname([pname fname], '' , '', TSTOOLdatapath, TSTOOLfilter);
	    sig=write(sig, newfilename);
	    cwd=pname;
	    set(loadhandle, 'UserData', {TSTOOLfilter cwd});
	    if ~isempty(newfilename)
	      datafiles=linsert(newfilename, handles,datafiles,sig,0);
	      if ploteverytime, tsplot(newfilename, plotmode); end
	    end
	  end
	 case 'impnetcdf'
	  [fname, pname] = uigetfile( fullfile(cwd, '*.nc'), 'Load NetCDF file', fileselboxXY(1), fileselboxXY(2));
	  if fname ~= 0
	    cwd=pname;
	    set(loadhandle, 'UserData', {TSTOOLfilter cwd});
	    sig = signal([pname fname], 'nc');
	    newfilename = newname([pname fname], '' , '', TSTOOLdatapath, TSTOOLfilter);
	    sig=write(sig, newfilename);
	    if ~isempty(newfilename)
	      datafiles=linsert(newfilename, handles,datafiles,sig,0);
	      if ploteverytime, tsplot(newfilename, plotmode); end
	    end
	  end
	  
	 case 'impmat'
	  [fname, pname] = uigetfile (fullfile(cwd, '*.mat'), ['Load' ...
		    ' Matlab-File'], fileselboxXY(1), ...
				      fileselboxXY(2));
	  if fname ~= 0
	    cwd=pname;
	    set(loadhandle, 'UserData', {TSTOOLfilter cwd});
	    var_cell = inputdlg ('Enter variable to load:', ...
				 'Loading data from matlab file');
	    varname = var_cell{1};
	    file = fullfile (pname, fname);
	    dat = load (file);
	    eval (['data = dat.' varname ';']);
	    
	    % now that we finally have extracted the data from our
            % mat file, we can construct our signal object
	    sig = signal(data);
	    newfilename = newname([pname fname], '', '', TSTOOLdatapath, ...
				  TSTOOLfilter);
	    sig = write (sig, newfilename);
	    if ~isempty(newfilename)
	      datafiles=linsert(newfilename, handles,datafiles,sig,0);
	      if ploteverytime, tsplot(newfilename, plotmode); end
	    end
	  end
	  
	 case 'embparm'
	  % Check for new embedding parameters
	  title = 'Reconstruction parameters';
	  lineNo = 1;
	  prompt  = {'Embedding dimension','Delay in samples', 'Shift in samples' ...
		    , 'Window type'};
	  defaults=get(recopthandle,'UserData');
		%	  defaults = { num2str(dim), num2str(delay), num2str(shift), recwindowtype};
	  answer  = sinputdlg(prompt,title,lineNo,defaults,{'','in units','in units',''});
	  if ~isempty(answer)
	    dim = str2num(answer{1});
	    delay = str2num(answer{2});
	    shift = str2num(answer{3});
	    recwindowtype = answer{4};
	    if (dim >= 0) & (delay >= 0)
	      set(recopthandle, 'UserData', answer(1:8));	
	    else
	      text = ['Embedding parameters not correct'];
	      warndlg(text);
	    end
	  end
	 case 'optmisc'
	  title = 'File and directory options';
	  lineNo = 1;
	  prompt  = {'Default file extension', 'Current working directory'};
	  defaults = get(loadhandle,'UserData');
	  answer  = inputdlg(prompt,title,lineNo,defaults);
	  if ~isempty(answer)
	    set(lallhandle, 'UserData', answer{1} );
	    TSTOOLfilter=answer{1};
	    if exist(answer{2}, 'dir') 		 % is a directory ?
	      cwd = answer{2};
	      set(loadhandle, 'UserData', answer);
	      set(fighandle, 'Name',['tstool - ' cwd]);
	    else
	      text = [answer{2}  ' seems not to be a directory'];
	      warndlg(text);
	    end
	  end
	 case 'togglepreviewmode'
	  nr = varargin{1};
	  if strcmp(get(plotmodehandle(nr), 'Checked'), 'on')
	    set(plotmodehandle(nr), 'Checked', 'off');
	  else
	    set(plotmodehandle(3-nr), 'Checked', 'off');	% switch other mode off
	    set(plotmodehandle(nr), 'Checked', 'on');
	  end
	 case 'sine'	
	  title = 'Sinus generator options';
	  lineNo = 1;
	  prompt  = { 'Length (in seconds)' , ...
		      'Frequency (Hz)', ...
		      'Amplitude', ...
		      'Samplerate (Hz)', ...
		      'Start Phase (rad)' };
	  %defaults = { 1, 1000, 1, 8000, 0 };
	  handle = findobj(fighandle, 'Tag', action);
	  defaults = get(handle, 'UserData');		
	  answer  = inputdlg(prompt,title,lineNo,defaults);
	  if ~isempty(answer)
	    len = str2num( answer{1} );
	    freq = str2num( answer{2} );
	    ampl = str2num( answer{3} );
	    rate = str2num( answer{4} );
	    phase = str2num( answer{5} );
	    set(handle, 'UserData', answer);
	    newfilename = newname(action, '' , '', TSTOOLdatapath, TSTOOLfilter);
	    sig = gen(action, len, freq, rate, ampl, phase);
	    sig=write(sig, newfilename);
	    if ~isempty(newfilename)
	      datafiles=linsert(newfilename, handles,datafiles,sig,0);
	      if ploteverytime, tsplot(newfilename, plotmode); end
	    end
	  end
	 case {'fnoise' , 'unoise', 'gnoise'}
	  title = 'Noise generator options';
	  lineNo = 1;
	  prompt  = {'Length (in seconds)', 'Samplerate (Hz)' };
	  handle = findobj(fighandle, 'Tag', action);
	  defaults = get(handle, 'UserData');
	  %defaults = { '1' ,  '8000'};
	  answer  = inputdlg(prompt,title,lineNo,defaults);
	  if ~isempty(answer) 
	    len = str2num( answer{1} );
	    rate = str2num( answer{2} ); 
	    set(handle, 'UserData', answer);
	    newfilename = newname(action, '' , '', TSTOOLdatapath, TSTOOLfilter);
	    sig = gen(action, len, 0, rate);
	    sig=write(sig, newfilename);
	    if ~isempty(newfilename)
	      datafiles=linsert(newfilename, handles,datafiles,sig,0);
	      if ploteverytime, tsplot(newfilename, plotmode); end
	    end
	  end
	 case 'dynsys'
	  title = 'Parameters for dynamical system integration';
	  lineNo = 1;
	  prompt  = {'Name of function to integrate' , ...
		     'Initial conditions', ...
		     'Length (in seconds)', ...
		     'Samplerate (Hz)', ...
		     'Solver'};
	  %defaults = {'lorenz', '[0; 0.01; -0.01]', 20, 100, 'ode45'};
	  handle = findobj(fighandle, 'Tag', action);
	  defaults = get(handle, 'UserData');
	  answer  = inputdlg(prompt,title,lineNo,defaults);
	  if ~isempty(answer)
	    fn = answer{1};
	    ic = str2num( answer{2} );
	    len = str2num( answer{3} );
	    rate = str2num( answer{4} );
	    solver = answer{5};
	    set(handle, 'UserData', answer);
	    newfilename = newname(fn, '' , '', TSTOOLdatapath, TSTOOLfilter);
	    sig = odesolver(fn, ic, len, rate, solver);
	    sig=write(sig, newfilename);
	    if ~isempty(newfilename)
	      datafiles=linsert(newfilename, handles,datafiles,sig,0);
	      if ploteverytime, tsplot(newfilename, plotmode); end
	    end
	  end
	 case 'genbyode'
	  [handle, figure] = gcbo;
	  defaults = get(handle, 'UserData');	
	  title = ['Generate ' defaults{1} ' data'];
	  lineNo = 1;
	  prompt  = {'Name of function to integrate' , ...
		     'System parameters', ...
		     'Initial conditions', ...
		     'Length (in seconds)', ...
		     'Samplerate (Hz)'};	
	  answer  = inputdlg(prompt,title,lineNo,defaults);
	  if ~isempty(answer)
	    system = answer{1};
	    params = str2num(answer{2});
	    ic = str2num( answer{3} );
	    len = str2num( answer{4} );
	    rate = str2num( answer{5} );
	    set(handle, 'UserData', answer);
	    newfilename = newname(system, '' , '', TSTOOLdatapath, TSTOOLfilter);
	    sig = genbyode(system, len, rate, ic, params); 
	    sig=write(sig, newfilename);
	    if ~isempty(newfilename)
	      datafiles=linsert(newfilename, handles,datafiles,sig,0);
	      if ploteverytime, tsplot(newfilename, plotmode); end
	    end
	  end
	 case 'genbymap'
	  [handle, figure] = gcbo;
	  defaults = get(handle, 'UserData');	
	  title = ['Generate ' defaults{1} ' data'];
	  lineNo = 1;
	  prompt  = {'Name of map' , ...
		     'Parameters', ...
		     'Initial conditions', ...
		     'Length (in samples)'};	
	  answer  = inputdlg(prompt,title,lineNo,defaults);
	  if ~isempty(answer)
	    system = answer{1};
	    params = str2num(answer{2});
	    ic = str2num( answer{3} );
	    len = str2num( answer{4} );
	    set(handle, 'UserData', answer);
	    newfilename = newname(system, '' , '', TSTOOLdatapath, TSTOOLfilter);
	    sig = genbymap(system, len, ic, params); 
	    sig=write(sig, newfilename);
	    if ~isempty(newfilename)
	      datafiles=linsert(newfilename, handles,datafiles,sig,0);
	      if ploteverytime, tsplot(newfilename, plotmode); end
	    end
	  end							
	 case 'editmacro' 
	  if exist(fullfile(TSTOOLdatapath, 'scripts', 'macro.m'))
	    l = readascii(fullfile(TSTOOLdatapath, 'scripts', 'macro.m'));
	    title = 'macro.m';
	    lineNo = 12;
	    prompt = {'Edit commands :'};
	    answer  = inputdialog(prompt,title,lineNo, {char(l)});
	    if ~isempty(answer)
	      l = list(cellstr(answer{1}));
	      writeascii(fullfile(TSTOOLdatapath, 'scripts', 'macro.m'), l);
	    end
	  else
	    warning('No macro defined yet');
	  end
	 case 'renamemacro'
	  if exist(fullfile(TSTOOLdatapath, 'scripts', 'macro.m'))
	    title = ['Rename macro'];
	    lineNo = 1;
	    prompt  = {'New filename '};
	    defaults = {'macro.m'};
	    answer  = inputdlg(prompt,title,lineNo,defaults);
	    if ~isempty(answer)
	      currentname = fullfile(TSTOOLdatapath, 'scripts', 'macro.m');
	      newfilename = fullfile(TSTOOLdatapath, 'scripts', answer{1});
	      status = copyfile(currentname, newfilename);
	      if status
		delete(currentname);
	      end
	    end
	  else
	    warning('No macro defined');
	  end		
	 case 'usage'
	  warndlg('Sorry, not yet supported');
	 otherwise
	  % Hier werden alle Callbacks behandelt, die einen Filenamen benötigen
	  currentfile =fullfile(TSTOOLdatapath,[getlastentry(datafiles,get(lboxhandle,'Value')) TSTOOLfilter]);
	  if  ~isempty(currentfile)
	    if exist(currentfile)
	      [path,name,ext,ver] = fileparts(currentfile);
	      switch action
	       case 'plot'
		tsplot(currentfile);
	       case 'info'
		desview(currentfile);
	       case {'spec', 'int', 'diff', 'db', 'abs', 'surrogate1', 'surrogate2', 'surrogate3', ...
		     'norm2', 'center', 'reverse', 'rang', 'histo', 'power', 'fft'}
		through(action, currentfile);		% Methoden ohne zusaetzliche Parameter, In->Out
	       case 'macro'
		 macfile=get(handles.macro,'UserData');
		 if exist(fullfile(TSTOOLdatapath, 'scripts', [macfile '.m']))
		  clear('macro'); 	% make shure newest macro.m is executed
		  through(macfile, currentfile);	
		else
		  warndlg('No macro defined yet');
		end
	      case 'save'
		
		[fname, pname] = uiputfile([getlastentry(datafiles,get(lboxhandle,'Value')) TSTOOLfilter], ...
		    'Save data file', fileselboxXY(1), ...
		    fileselboxXY(2));
		
		if fname ~= 0
		  fullname = fullfile(pname, fname);
		  s=signal(currentfile);
		  write(s,fullname);
		end
		
	       case 'embed'
		title = 'Reconstruction parameters';
		lineNo = 1;
		
		prompt={'Embedding dimension', ...
			'Delay in samples', ...
			'Shift in samples', ...
			'Window type'};
		defaults=get(recopthandle,'UserData');
%		defaults = { num2str(dim), num2str(delay), num2str(shift), recwindowtype,0,0,0,0};
		answer  = sinputdlg(prompt,title,lineNo,defaults,{'','in units','in units',''});
		if ~isempty(answer)
		  dim = str2num(answer{1});
		  delay = str2num(answer{2});
		  shift = str2num(answer{3});
		  recwindowtype = answer{4};
		  if (dim >= 0) & (delay >= 0)
		    set(recopthandle, 'UserData', answer(1:8));	
		  else
		    text = ['Embedding parameters not correct'];
		    warndlg(text);
		  end
		  through('embed', currentfile, dim, delay, shift, recwindowtype,answer{5:8},'unitsamps');
		end
		
	       case 'mixembed'
		title    = 'Reconstruction parameters';
		lineNo   = 1;
		prompt   = {'Embedding dimension 1', ...
			    'Embedding dimension 2', ...
			    'Delay 1 in samples', ...
			    'Delay 2 in samples', ...
			    'Shift in samples', ...
			    'Window type'};
		
		defaults = {'1', '1', '1', '1', '0', 'Rect'};
		answer   = inputdlg (prompt, title, lineNo, defaults);
		
		if ~isempty(answer)
		  dim1   = str2num(answer{1});
		  dim2   = str2num(answer{2});
		  delay1 = str2num(answer{3});
		  delay2 = str2num(answer{4});
		  shift  = str2num(answer{5});
		  recwindowtype = answer{6};
		  if ((dim1 >= 0) & ...
		      (dim2 >= 0) & ...
		      (delay1 >= 1) & ...
		      (delay2 >= 1) & ...
		      (shift >= 0))
		    set (recopthandle, 'UserData', answer);
		  else
		    text = ['Embedding parameters not correct'];
		    warndlg (text);
		  end
		  through ('mixemb', currentfile, dim1, dim2, ...
			   delay1, delay2, shift);
		end
		
	       case    {'trev','tc3'}
		title='Surrogate data test parameters';
		prompt = {'Delay (in Samples)','Number of tests', ...
			  'Surrogate generation method (1,2,3)'};
		defaults=get(handles.trevtc3,'UserData');
%		defaults={'10','1','1',0,0,0};
		answer=sinputdlg(prompt,title,1,defaults,{'in units','',''});
		if ~isempty(answer)
		  set(handles.trevtc3,'UserData',answer);
		  through(action,currentfile,str2num(answer{1}),str2num(answer{2}),str2num(answer{3}), ...
		      answer{4:end},'unitsamps');
		end 
	       case 'surro'     
		title='Surrogate data test parameters';
		prompt = {'Number of tests', ...
			  'Surrogate generation method (1,2,3)', ...
			  'Function'};
		defaults=get(handles.surrogate,'UserData');
%		defaults={'10','1', ...
%			  'corrsum(embed(s,3,0,0), -1, 0.1, 20, 32, 0);'};
		
		answer  = inputdlg(prompt,title,1,defaults);		
		if ~isempty(answer)
		  set(handles.surrogate,'UserData',answer);
		  through('surrogate_test',currentfile,str2num(answer{1}), ...
			  str2num(answer{2}),answer{3});
		end
		
	       case 'stts'
		title = 'Reconstruction parameters';
		prompt  = {'Number of spatial neighbours','Number of temporal neighbours (in the past)', ...
			   'Spatial shift (= spatial delay)', 'Temporal delay'};
		defaults = { '1' , '1', '1', '1' };
		answer  = inputdlg(prompt,title,1,defaults);
		if ~isempty(answer)
		  I= str2num( answer{1} );J= str2num( answer{2} );
		  K= str2num( answer{3} );L= str2num( answer{4} );
		  through('stts', currentfile, I,J,K,L);
		end
	       case 'multires'
		title = 'Multiresolution analysis';
		lineNo = 1;
		prompt  = {'Depth'};
		defaults = {'3' };
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  depth = str2num(answer{1});
		  through('multires', currentfile, depth);
		end						
	       case 'spec2'
		title = 'Spectrogram';
		lineNo = 1;
		prompt  = {'Window length (in samples)','Window type', 'Window shift (in samples)'};
		defaults = {'512', windowtype, '128', 0, 0, 0};
		answer  = sinputdlg(prompt,title,lineNo,defaults,{'in units','','in units'});
		if ~isempty(answer)
		  flen= str2num( answer{1} );
		  ftype = answer{2}
		  vorschub = str2num( answer{3} );
		  through('spec2', currentfile, flen, ftype, vorschub,answer{4:end},'unitsamps');
		end
	       case 'scalogram'
		title = 'Scalogram';
		lineNo = 1;
		prompt  = {'Start scale', 'End scale', 'Scalestep', 'Tlength'};
		defaults = {'1', '5', '0.5', '8'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  smin = str2num(answer{1});
		  smax = str2num(answer{2});
		  sstep = str2num(answer{3});
		  tlen = str2num(answer{4});
		  through('scalogram', currentfile,smin, smax, sstep, tlen);
		end	
	       case 'aok'
		title = 'AOK TFR';
		lineNo = 1;
		prompt  = {'Window length (in samples)','FFT Length', 'Window shift', 'Kernel Volume'};
		defaults = { '32', '256', '1' '2'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  wlen= str2num(answer{1});
		  flen = str2num(answer{2});
		  vorschub = str2num( answer{3});
		  volume = str2num( answer{4});
		  through('aok', currentfile, wlen, flen, vorschub, volume);
		end
		
	       case 'nnstatistic'
		title = 'Nearest neighbor statistic';
		lineNo = 1;
		prompt  = {'Number of nearest neighbours', 'Exclude next and last samples'};
		defaults = { '8', '10' };
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  nnr = str2num( answer{1} );
		  exclude = str2num(answer{2});
		  through2('nearneigh', currentfile, nnr, exclude);
		end
	       case 'ret_time'
		title = 'Return times';
		lineNo = 1;
		prompt  = {'Number of nearest neighbours', 'Maximal return time', 'Exclude next and last samples'};
		defaults = { '8', '128', '10' };
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  nnr = str2num( answer{1} );
		  maxT = str2num( answer{2} );
		  exclude = str2num(answer{3});
		  through2('return_time', currentfile, nnr, maxT, exclude);
		end
	       case 'density'
		title = 'Reciprocal local density';
		lineNo = 1;
		prompt  = {'Number of nearest neighbours', 'Exclude next and last samples'};
		defaults = { '8', '10' };
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  nnr = str2num( answer{1} );
		  exclude = str2num(answer{2});
		  through2('localdensity', currentfile, nnr, exclude);
		end
	       case 'largelyap'
		title = 'Estimate largest lyapunov exponent';
		lineNo = 1;
		prompt  = {'Number of reference points', 'Length of prediction', 'Number of nearest neighbours', 'Exclude next and last samples'};
		defaults = { '100', '20', '1', '10' };
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  n = str2num( answer{1} );
		  len = str2num( answer{2} );
		  nnr = str2num( answer{3} );
		  exclude = str2num(answer{4});
		  through2('largelyap', currentfile, n, len, exclude, nnr);
		end
	       case 'corrsum'
		title = 'Estimate correlation dimension';
		lineNo = 1;
		prompt  = {'Number of reference points', 'Relative search radius', 'Exclude next and last samples', 'Number of bins'};
		defaults = { '1000', '0.1', '20', '32' };
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  n = str2num( answer{1} );
		  range = str2num( answer{2} );
		  exclude = str2num(answer{3});
		  bins = str2num( answer{4} );
		  through2('corrsum', currentfile, n, range, exclude, bins);
		end								
	       case 'corrsum2'
		title = 'Estimate correlation dimension';
		lineNo = 1;
		prompt  = {'Number of pairs', 'Relative search radius', 'Exclude next and last samples', 'Number of bins'};
		defaults = { '[1000 100 2000]', '0.05', '20', '32' };
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  npairs = str2num( answer{1} );
		  range = str2num( answer{2} );
		  exclude = str2num(answer{3});
		  bins = str2num( answer{4} );
		  through2('corrsum2', currentfile, npairs, range, exclude, bins);
		end
		
	       case  'tak_est'

		title = 'Takens estimator for corr. dim';
		lineNo = 1;
		prompt  = {'Number of reference points', 'Relative search radius', 'Exclude next and last samples'};
		defaults = { '1000', '0.05', '20'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  n = str2num( answer{1} );
		  range = str2num( answer{2} );
		  exclude = str2num(answer{3});
		  sig = signal(currentfile);
		  est = takens_estimator(sig, n, range, exclude);
		  msgbox({'Takens estimator of corr. dimension' ['D2 = ' num2str(est)]}, name);
		end
	       case 'infodim2'
		title = 'Estimate information dimension';
		lineNo = 1;
		prompt  = {'Number of reference points', 'Maximal numbers of neighbors', 'Exclude next and last samples'};
		defaults = { '1000', '128', '20'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  n = str2num( answer{1} );
		  kmax = str2num( answer{2} );
		  exclude = str2num(answer{3});						
		  through2('infodim2', currentfile, n, kmax, exclude);
		end
	       case 'fracdims'
		title = 'Estimate fractal dimension spectrum';
		lineNo = 1;
		prompt  = {'Number of reference points', 'Minimal numbers of neighbors', ...
			   'Maximal numbers of neighbors', 'Starting moment order', 'End moment order', ...
			   'Exclude next and last samples'};
		defaults = { '1000', '8', '128', '-4' , '4', '0'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  nref = str2num( answer{1} );
		  kmin = str2num( answer{2} );
		  kmax = str2num( answer{3} );
		  gmin = str2num( answer{4} );
		  gmax = str2num( answer{5} );
		  exclude = str2num(answer{6});						
		  through2('fracdims', currentfile, kmin, kmax, nref, gmin, gmax, exclude);
		end												   
	       case 'boxdim'
		title = 'Estimate boxcounting dimension';
		lineNo = 1;
		prompt  = {'Maximal number of bins'};
		defaults = {'128'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  bins = str2num( answer{1} );
		  through('boxdim', currentfile, bins);
		end
	       case 'infodim'
		title = 'Estimate information dimension';
		lineNo = 1;
		prompt  = {'Maximal number of bins'};
		defaults = {'128'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  bins = str2num( answer{1} );
		  through('infodim', currentfile, bins);
		end   
	       case 'corrdim'
		title = 'Estimate correlation dimension';
		lineNo = 1;
		prompt  = {'Maximal number of bins'};
		defaults = {'128'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  bins = str2num( answer{1} );
		  through('corrdim', currentfile, bins);
		end     
	       case 'localconstpred'
		title = 'Local constant prediction';
		lineNo = 1;
		prompt  = {'Length of prediction', 'Number of nearest neighbours', 'Stepsize', 'Averaging mode'};
		defaults = {'100', '1', '1', '0'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  len = str2num(answer{1});
		  nnr = str2num(answer{2});
		  step = str2num(answer{3});
		  mode = str2num(answer{4});
		  through('predict', currentfile, len, nnr, step, mode);
		end
	       case 'scalarpredict'
		title = 'Local constant prediction';
		lineNo = 1;
		prompt  = {'Dimension', 'Delay', 'Length of prediction', 'Number of nearest neighbours', 'Averaging mode'};
		defaults = {num2str(dim), num2str(delay), '100', '1', '0'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  dim = str2num(answer{1});
		  delay = str2num(answer{2});
		  len = str2num(answer{3});
		  nnr = str2num(answer{4});
		  mode = str2num(answer{5});
		  through('predict2', currentfile, dim, delay, len, nnr, mode);
		end						
	       case 'acf'
		title = 'Autocorrelation';
		prompt  = {'FFT Length (in samples)'};
		defaults = { '256',0};
		answer  = sinputdlg(prompt,title,1,defaults,{'in units'});
		if ~isempty(answer)
		  len= str2num( answer{1} );
		  through('acf', currentfile, len,answer{2},'unitsamps');
		end
	       case 'ccf'
		title = 'Cross-correlation';
		prompt  = {'Length'};
		defaults = { '256',0};
		answer  = sinputdlg(prompt,title,1,defaults,{'in units'});
		secfilename = secfile(get(handles.lboxhandle, 'String'));
		if ~isempty(secfilename)
		  len= str2num( answer{1} );
		  through3('ccf', currentfile, secfilename,len);
		end
	       case 'gmi'
		title = 'Generalized mutual information';
		lineNo = 1;
		prompt  = {'Length', 'Number of reference points', 'Correlation radius', 'Maximal number of nearest neighbours'};
		defaults = { '32', '500', '0.05', '512' };
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  len = str2num( answer{1} );
		  n = str2num( answer{2} );
		  range = str2num( answer{3} );
		  nnr = str2num( answer{4} );
		  through('gmi', currentfile, 2, range, nnr, len, n);
		end
	       case 'cao'
		title = 'Estimate minimum embedding dimension';
		lineNo = 1;
		prompt  = {'Maxdim', 'Delay', 'Maximal number of nearest neighbours', 'Number of reference points'};
		defaults = { '10', '1', '1', '100' };
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  maxdim = str2num( answer{1} );
		  tau = str2num( answer{2} );
		  nnr = str2num( answer{3} );
		  nref = str2num( answer{4} );
		  through('cao', currentfile, maxdim, tau, nnr, nref);
		end
	       case  'polmodel'
		title = 'Polynom selection';
		lineNo = 1;
		prompt  = {'Maximal degree', 'Fraction of test/training samples'};
		defaults = { '5', '0.5'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  maxgrad = str2num( answer{1} );
		  fracref = str2num( answer{2} );
		  sig = signal(currentfile);
		  dat = data(sig);
		  pol = pauswahl(dat(1:end-1,:), dat(2:end,end), fracref, maxgrad);
		  %msgbox
		end							
	       case 'poincare'
		title = 'Poincare section';
		prompt  = {'Index of reference sample (2..END-1)'};
		defaults = { '2'};
		answer  = inputdlg(prompt,title,1,defaults);
		if ~isempty(answer)
		  ref = str2num( answer{1} );
		  through('poincare', currentfile, ref);
		end
	       case 'amutual'
		title = 'Auto mutual information';
		lineNo = 1;
		prompt  = ...
		    {'Length in samples'};
		defaults = { '128',0};
		answer  = sinputdlg(prompt,title,lineNo,defaults,{'in units'});
		if ~isempty(answer)
		  len= str2num( answer{1} );
		  through('amutual', currentfile, len,answer{2},'unitsamps');
		end
	       case 'cubicsplineinterp'
		title = 'Cubic spline interpolation';
		lineNo = 1;
		prompt  = {'Upsampling factor'};
		defaults = {'2'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  factor = str2num( answer{1} );
		  through('upsample', currentfile, factor, 'spline');
		end
	       case 'akimasplineinterp'
		title = 'Akima spline interpolation';
		lineNo = 1;
		prompt  = {'Upsampling factor'};
		defaults = {'2'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  factor = str2num( answer{1} );
		  through('upsample', currentfile, factor, 'akima');
		end
	       case 'fftinterp'
		title = 'Increase sample rate';
		lineNo = 1;
		prompt  = {'Upsampling factor'};
		defaults = {'2'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  factor = str2num( answer{1} );
		  through('upsample', currentfile, factor , 'fft');
		end
	       case 'movav'
		title = 'Moving average';
		lineNo = 1;
		prompt  = {'Window length (in samples)'};
		defaults = { '5',0};
		answer  = sinputdlg(prompt,title,lineNo,defaults,{'in units'});
		if ~isempty(answer)
		  len= str2num( answer{1} );
		  through('movav', currentfile, len,answer{2},'unitsamps');
		end
	       case 'medianfilt'
		title = 'Moving median filter';
		lineNo = 1;
		prompt  = {'Window length (in samples)'};
		defaults = { '5',0};
		answer  = sinputdlg(prompt,title,lineNo,defaults,{'in units'});
		if ~isempty(answer)
		  len= str2num( answer{1} );
		  through('medianfilt', currentfile, len,answer{2},'unitsamps');
		end
	       case 'level_adaption'
		title = 'Automatic level adaption';
		lineNo = 1;
		prompt  = {'Time constants', 'Dynamic limit ','Threshold'};
		defaults = {'[500 300 100]' , '10', '0.0005'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)
		  time_constants = str2num(answer{1});
		  limit = str2num(answer{2});
		  thresh = str2num(answer{3});
		  through('level_adaption', currentfile, time_constants, limit, thresh);
		end						
		%                         case 'plosivity'
		% 						title = 'Plosivity';
		% 						lineNo = 1;
		% 						prompt  = {'Backward length (in samples)', 'Forward length (in samples)','Threshold'};
		% 						defaults = {'5' , '2', '0'};
		% 						answer  = inputdlg(prompt,title,lineNo,defaults);
		% 						if ~isempty(answer)
		% 							blen = str2num( answer{1} );
		%                                 flen = str2num( answer{2} );
		%                                 thresh = str2num( answer{3} );
		% 							through('plosivity', currentfile, blen, flen, thresh);
		%                         end
	       case 'histo'
		title = 'Histogram';
		prompt  = {'Number of bins'};
		defaults = { '100'};
		answer  = inputdlg(prompt,title,1,defaults);
		if ~isempty(answer)
		  bins = str2num( answer{1} );
		  through('histo', currentfile, bins);
		end
	       case 'trend'
		title = 'Trend correction';
		lineNo = 1;
		prompt  = {'Average length (in samples)'};
		defaults = {'9',0};
		answer  = sinputdlg(prompt,title,lineNo,defaults,{'in units'});
		if ~isempty(answer)
		  len= str2num(answer{1} );
		  through('trend', currentfile, len,answer{2},'unitsamps');
		end
	       case 'max'
		[y, yu, xpos, xu] = max(signal(currentfile));
		msgbox({'Maximum of signal :' ['y = ' num2str(y) ' ' label(yu)] ['x = ' num2str(xpos) ' ' label(xu)] }, name);
	       case 'min'
		[y, yu, xpos, xu] = min(signal(currentfile));
		msgbox({'Minimum of signal' ['y = ' num2str(y) ' ' label(yu)] ['x = ' num2str(xpos) ' ' label(xu)] }, name);		
	       case 'firstmin'
		[xpos, u] = firstmin(signal(currentfile));
		if ~isempty(xpos)
		  msgbox({'First local minimum at ' ['x = ' num2str(xpos) ' ' label(u)]}, name);		
		else
		  msgbox('No local minimum found', name);
		end						
	       case 'firstmax'
		[xpos, u] = firstmax(signal(currentfile));
		if ~isempty(xpos)
		  msgbox({'First local maximum at ' ['x = ' num2str(xpos) ' ' label(u)]}, name);		
		else
		  msgbox('No local maximum found', name);
		end	
	       case 'firstzero'
		[xpos, u] = firstzero(signal(currentfile));
		if ~isempty(xpos)
		  msgbox({'First zero crossings at ' ['x = ' num2str(xpos) ' ' label(u)]}, name);		
		else
		  msgbox('No zero crossings found', name);
		end	
	       case 'rms'
		sig = signal(currentfile);
		u = yunit(sig); v = data(sig);
		v = sqrt(sum(v.*v)/dlens(sig,1));
		msgbox({'Root mean square of signal' ['RMS = ' num2str(v) ' ' label(u)]}, name);
	       case 'mean'
		sig = signal(currentfile);
		u = yunit(sig); v = mean(data(sig));
		msgbox({'Mean of signal' [num2str(v) ' ' label(u)]}, name);
	       case 'std'
		sig = signal(currentfile);
		u = yunit(sig); v = std(data(sig));
		msgbox({'Standard deviation of signal' [num2str(v) ' ' label(u)]}, name);
	       case 'swap'
		sig = signal(currentfile);	
		if ~isempty(sig)
		  switch ndim(sig)
		   case 1
		    warndlg('Signal has only one dimension, swapping makes no sense');
		    return
		   case 2
		    sig = swap(sig);
		   otherwise
		    title = ['Swap dimensions'];
		    lineNo = 1;
		    prompt  = {'Swap dimension Nr. ' , 'with dimension Nr. '};
		    defaults = { '1' , '2' };
		    answer  = inputdlg(prompt,title,lineNo,defaults);
		    if ~isempty(answer)
		      first = str2num(answer{1});
		      second = str2num(answer{2});
		      if (first<ndim(sig)) | (second<ndim(sig))
			sig = swap(sig, dim1, dim2);
		      else
			warndlg('Input exceed signal''s dimensions');
			return
		      end
		    else
		      return
		    end
		  end
		  if ~isempty(sig)
		    newfilename = newname(currentfile, 'swap_', '', TSTOOLdatapath);
		    sig=write(sig, newfilename);
		    if ~isempty(newfilename)
		      datafiles=linsert(newfilename, handles,datafiles,sig,0);
		      if ploteverytime
%			if strcmp(plotmode, 'small')
%			  status = view(sig);
%			  zoom off
%			  rotate3d off
%			else
			  tsplot(newfilename, plotmode)
%			end
		      end
		    end
		  else
		    return;
		  end
		end
	       case 'norm1'
		title = ['Fit signal values to interval'];
		lineNo = 1;
		prompt = { 'Lower limit', 'Upper limit'};
		defaults  = {'0','1'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)	
		  low = str2num(answer{1});
		  upp = str2num(answer{2});
		  if upp > low
		    through('norm1', currentfile, low, upp);
		  else
		    warndlg('Wrong interval boundaries given');
		  end
		end					
	       case 'scale'
		title = ['Scale signal'];
		lineNo = 1;
		prompt = { 'Scaling Factor'};
		defaults  = {'1'};
		answer  = inputdlg(prompt,title,lineNo,defaults);
		if ~isempty(answer)	
		  fct = str2num(answer{1});
		  if fct ~= 1
		    through('scale', currentfile, fct);
		  else
		    warndlg('Nothing to do');
		  end
		end
	       case 'cut'
		sig = signal(currentfile);	
		if ~isempty(sig)
		  title = ['Cut signal'];
		  lineNo = 1;
		  prompt  = {'Along dimension' , 'from start index' , 'to stop index'};
		  defaults = { '1' , '1', num2str(dlens(sig, 1)) };
		  answer  = inputdlg(prompt,title,lineNo,defaults);
		  if ~isempty(answer)
		    dim = str2num(answer{1});
		    start = str2num(answer{2});
		    stop  = str2num(answer{3});
		    if (dim > ndim(sig)) | (start > dlens(sig,dim)) | (stop<start)
		      warndlg('Wrong parameters given');
		      return
		    end
		    sig = cut(sig, dim, start, stop);
		    newfilename = newname(currentfile, 'cut_', '', TSTOOLdatapath);
		    sig=write(sig, newfilename);
		    if ~isempty(newfilename)
		      datafiles=linsert(newfilename, handles,datafiles,sig,0);		      
		      if ploteverytime
%			if strcmp(plotmode, 'small')
%			  status = view(sig);
%			  zoom off
%			  rotate3d off
%			else
			  tsplot(newfilename, plotmode)
%			end
		      end
		    end
		  end
		end
	       case 'split'
		sig = signal(currentfile);
		if ~isempty(sig)
		  if ndim(sig)~=2
		    warndlg([name ' is not a multichannel signal']);
		    return
		  end
		  for n=1:dlens(sig,2)
		    sig2 = cut(sig, 2, n, n);
		    newfilename = newname(currentfile, '', num2str(n), TSTOOLdatapath);
		    if n==1
		      [dummy,fname,dummy,dummy]=fileparts(newfilename);
		    end
		    
		    sig2=write(sig2, newfilename);
		    if ~isempty(newfilename)
		      datafiles=linsert(newfilename, handles,datafiles,sig2,1);
		    end
		  end
		  datafiles=sortdatafiles(datafiles);
		  filllistbox(datafiles,handles,fname);
		  
		end
	       case 'pca'
		sig = signal(currentfile);
		if ndim(sig) ~= 2
		  warndlg('pca needs a signal with two dimensions as input')
		  return
		end
		if ~isempty(sig)
		  title = ['Principal component analysis'];
		  lineNo = 1;
		  prompt = { 'Preprocessing (norm, mean, raw)' , 'Cumulative Percentage' };
		  defaults  = {'norm' , '95'};
		  answer  = inputdlg(prompt,title,lineNo,defaults);
		  if ~isempty(answer)
		    mode = answer{1};
		    maxpercent= str2num(answer{2});
		    [rs, eigvals, eigvecs] = pca(sig, mode, maxpercent);	
		  else
		    return
		  end						
		  newfilename = newname(currentfile, ['pca_'], '', TSTOOLdatapath);
		  rs=write(rs, newfilename);
		  if ~isempty(newfilename)
		    datafiles=linsert(newfilename, handles,datafiles,rs,1);
		  end
		  newfilename = newname(currentfile, ['eigvals_'], '', TSTOOLdatapath);
		  eigvals=write(eigvals, newfilename);
		  if ~isempty(newfilename)
		    datafiles=linsert(newfilename, handles,datafiles,eigvals,1);
		  end
		  newfilename = newname(currentfile, ['eigvecs_'], '', TSTOOLdatapath);
		  eigvecs=write(eigvecs, newfilename);
		  if ~isempty(newfilename)
		    datafiles=linsert(newfilename, handles,datafiles,eigvecs,0);
		  end
		end
	       case 'archetypes'
		sig = signal(currentfile);
		if ndim(sig) ~= 2
		  warndlg('arch needs a signal with two dimensions as input')
		  return
		end
		if ~isempty(sig)
		  title = ['Archetypal analysis'];
		  lineNo = 1;
		  prompt = {'Number of archetypes', 'Preprocessing (norm, mean, raw)'};
		  defaults  = {'5' , 'norm'};
		  answer  = inputdlg(prompt,title,lineNo,defaults);
		  if ~isempty(answer)
		    na = answer{1};
		    mode = answer{2};
		    [rs, archetypes] = arch(sig, na, mode);	
		  else
		    return
		  end						
		  newfilename = newname(currentfile, ['arch_'], '', TSTOOLdatapath);
		  rs=write(rs, newfilename);
		  if ~isempty(newfilename)
		    datafiles=linsert(newfilename, handles,datafiles,rs);
		  end
		  newfilename = newname(currentfile, ['archetypes_'], '', TSTOOLdatapath);
		  archetypes_=write(archetypes_, newfilename);
		  if ~isempty(newfilename)
		    datafiles=linsert(newfilename, handles,datafiles,archetypes_);
		  end
		end						
	       case {'plus', 'minus'}
		secfilename = secfile(get(handles.lboxhandle, 'String'));
		if ~isempty(secfilename)
		  through3(action, currentfile, secfilename);
		end
	       case 'merge'
		secfilename = secfile(get(handles.lboxhandle, 'String'));
		if ~isempty(secfilename)
		  title = ['Merge two signals'];
		  lineNo = 1;
		  prompt  = {'Desired energy ratio (dB) :'};
		  defaults = {'0'};
		  answer  = inputdlg(prompt,title,lineNo,defaults);
		  if ~isempty(answer)
		    through3(action, currentfile, secfilename, str2num(answer{1}));
		  end
		end
	       case 'compare'
		secfilename = secfile(get(handles.lboxhandle, 'String'));
		if ~isempty(secfilename)
		  sig = signal(currentfile);
		  sig2 = signal(secfilename);
		  [path,name1,ext,ver] = fileparts(currentfile);
		  [path,name2,ext,ver] = fileparts(secfilename);
		  if compare(sig,sig2)
		    msgbox([name1 ' and ' name2 ' have equal values'],'Compare signals')
		  else
		    msgbox([name1 ' and ' name2 ' are different'],'Compare signals')
		  end
		end
		
		% Export to other formats	
	       case 'ascii'			
		[fname, pname] = uiputfile( [name '.dat'], 'Export file as ASCII file', fileselboxXY(1), fileselboxXY(2));
		if fname ~= 0
		  sig = signal(currentfile);
		  write(sig, [pname fname], 'ASCII');
		end
	       case 'mat'
		[fname, pname] = uiputfile( [name '.mat'], 'Export file as .mat file', fileselboxXY(1), fileselboxXY(2));
		if fname ~= 0
		  sig = signal(currentfile);
		  write(sig, [pname fname], 'MAT');
		end
	       case 'sipp'
		[fname, pname] = uiputfile( [name '.si'], 'Export file as si++ file', fileselboxXY(1), fileselboxXY(2));
		if fname ~= 0
		  sig = signal(currentfile);
		  write(sig, [pname fname], 'SIPP');
		end
	       case 'nld'
		[fname, pname] = uiputfile( [name '.nld'], 'Export file as .nld file', fileselboxXY(1), fileselboxXY(2));
		if fname ~= 0
		  sig = signal(currentfile);
		  write(sig, [pname fname], 'NLD');
		end
	       case 'da'
		sig = signal(currentfile);
		if ndim(sig) == 1
		  a = getaxis(sig);
		  %if delta(a) > 0								
		  %rate = 	1/delta(a);
		  %else
		  rate = 8000;	% use a default sampling rate
				%	end
				soundsc(data(sig), rate);
		else
		  error('Sound output only for scalar signals');
		end
	       case 'edit'
		hded1(currentfile);
	       case 'editplothint'	
		sig = signal(currentfile);
		sig = setplothint(sig, varargin{1});
		sig = addcommandlines(sig, 's = setplothint(s', varargin{1});
		write(sig, currentfile);
	       case 'editaxes'
		hded2(currentfile);
	       case 'newcomment'
		text = 'Newly added comment text';	
		[path,name,ext,ver] = fileparts(currentfile);
		sig = signal(currentfile);
		title = ['Comment for ' name];
		lineNo = 12;
		prompt = {'New comment text :'};
		answer  = inputdialog(prompt,title,lineNo, {comment(sig)});
		if ~isempty(answer)
		  sig = newcomment(sig, answer{1});
		  write(sig, currentfile);
		end
	       case 'makescript'
		title = ['Create script from commandlines'];
		lineNo = 1;
		prompt = { 'Script name'};
		defaults = {fullfile(TSTOOLdatapath, 'scripts', 'newscript.m')};
		answer  = inputdlg(prompt,title,lineNo, defaults);
		if ~isempty(answer)	
		  sig = signal(currentfile);	
		  makescript(sig, answer{1});                
		end
	       case 'makemacro'		% Makro erzeugen
		sig = signal(currentfile);
		makescript(sig, fullfile(TSTOOLdatapath, 'scripts', 'macro.m'));
	       case 'history'
		text = 'History';	
		[path,name,ext,ver] = fileparts(currentfile);
		sig = signal(currentfile);
		title = ['History for ' name];
		lineNo = 12;
		prompt = {''};
		answer  = inputdialog(prompt,title,lineNo, {history(sig)});
	       case 'Rm!!!'	% entfernt File aus der Liste und physikalisch !!!
		answer = questdlg(['Physically delete file ' currentfile], ...
				  'tstool','No');
		if  strcmp(answer, 'Yes')
		  delete(currentfile);
		  lremove(handles);
		  cla
		end
	       otherwise
		error(['There''s no entry in tscallback.m for a callback named : ' action]);
	      end 	% switch
	    else
	      warndlg(['File ' currentfile ' does not exist. Try rescan.']);
	      
	    end 	% if exist
	  else
	    warndlg('No file selected');
	  end 	% if not empty	
end
	catch 			
  		errordlg(lasterr, 'TSTOOL-Error', 'on');   	
end

% Ende der Hauptroutine tscallback

% Folgende Routine wird jedesmal aufgerufen, wenn durch die Verarbeitung eines Inputfiles
% (meistens der currentfile) ein neues Signal entsteht und dieses automatisch abgespeichert
% und dargestellt werden soll (entweder im Preview-Fenster (plotmode = 'small') oder in einem 
% eigenen Window (plotmode = 'large'))
%
% Zuerst wird das Inputfile geladen (die Routine bekommt nur den Namen des Inputfiles uebergeben)
% Falls alles stimmt, wird die gewuenschte Verarbeitung 

function [status, newfilename] = through(func, filename, varargin)	

global TSTOOLdatapath;

fighandle = findobj('Tag', 'TSTOOL');  % Filename ist in den UserData der TSTOOL-Figure gespeichert
handles = get(fighandle, 'UserData');	% structure with handles to the different menus/controls of the figure
answer=get(handles.loadhandle,'UserData');
TSTOOLfilter=answer{1};
cwd=answer{2};
datafiles=get(handles.lboxhandle,'UserData');
ploteverytime = 0;
if strcmp(get(handles.plotmodehandle(1), 'Checked'), 'on')
  plotmode = 'small';
  ploteverytime = 1;
end
	
if strcmp(get(handles.plotmodehandle(2), 'Checked'), 'on')
  plotmode = 'large';
  ploteverytime = 1;
end

newfilename = '';
sig = signal(filename);

if ~isempty(sig)
  if nargin>2 & strcmp(varargin{end},'unitsamps')
    varargin=varargin(1:end-1);
    for i=1:length(varargin)/2
      if varargin{i+length(varargin)/2}==1
	ax=getaxis(sig,varargin{i+length(varargin)/2});
	varargin{i}=round(varargin{i}/delta(ax));
      end    
    end
    newsig = feval(func, sig, varargin{1:end/2});
  else
    newsig = feval(func, sig, varargin{:});
  end
  if ~isempty(newsig)
    newfilename = newname(filename, [func '_'], '', TSTOOLdatapath);
    newsig=write(newsig, newfilename);
    if ~isempty(newfilename)
      [status,datafiles]=siglevel(newsig,datafiles);
      datafiles=sortdatafiles(datafiles);
      %			disp(['Dateiname >' name(newsig) '<']);
      filllistbox(datafiles,handles,name(newsig));
      %		        datafiles=linsert(newfilename, handles,datafiles);
      if ploteverytime
%	if strcmp(plotmode, 'small')
%	  view(newsig, 6);
%	  %set(gca, 'Visible', 'off');
%	  zoom off
%	  rotate3d off
%	else
	  tsplot(newfilename, plotmode)
%	end
      end
    end
  else
    return;
  end
else
  return;
end

% through2 ist ein Klone der Funktion through2, der dann eingesetzt wird, wenn
% das Eingangssignal veraendert mit ausgegenen wird, was z.B. bei den NN-basierten
% Routinen der Fall sein kann, wenn die NN-Vorverarbeitung mit im Signal abgelegt wird.

function [status, newfilename] = through2(func, filename, varargin)	

global TSTOOLdatapath;


fighandle = findobj('Tag', 'TSTOOL');  % Filename ist in den UserData der TSTOOL-Figure gespeichert
handles = get(fighandle, 'UserData');	% structure with handles to the different menus/controls of the figure
answer=get(handles.loadhandle,'UserData');
TSTOOLfilter=answer{1};
cwd=answer{2};
datafiles=get(handles.lboxhandle,'UserData');
ploteverytime = 0;
if strcmp(get(handles.plotmodehandle(1), 'Checked'), 'on')
  plotmode = 'small';
  ploteverytime = 1;
end

if strcmp(get(handles.plotmodehandle(2), 'Checked'), 'on')
  plotmode = 'large';
  ploteverytime = 1;
end

newfilename = '';
sig = signal(filename);

if ~isempty(sig)
  if isempty(optparams(sig))
    no_preprocessing_found = 1;
  else
    no_preprocessing_found = 0;
  end
  
  [newsig, sig] = feval(func, sig, varargin{:});
  
  if (~isempty(optparams(sig))) & no_preprocessing_found 
    sig=write(sig, filename);
  end
  
  if ~isempty(newsig)
    newfilename = newname(filename, [func '_'], '', TSTOOLdatapath);
    newsig=write(newsig, newfilename);
    if ~isempty(newfilename)
      [status,datafiles]=siglevel(newsig,datafiles);
      datafiles=sortdatafiles(datafiles);
      filllistbox(datafiles,handles,name(newsig));
      if ploteverytime
%	if strcmp(plotmode, 'small')
%	  view(newsig, 6);
%	  %set(gca, 'Visible', 'off');
%	  zoom off
%	  rotate3d off
%	else
	  tsplot(newfilename, plotmode)
%	end
      end
    end
  else
    return;
  end
else
  return;
end


% through3 ist wiederum ein Klone von through2, der dann verwendet
% wird, wenn das zweite Argument ein Filename eines zweiten Signals
% ist (fuer plus,minus und merge). (DE)

function [status, newfilename] = through3(func, filename, secfile, varargin)	

global TSTOOLdatapath;

fighandle = findobj('Tag', 'TSTOOL');  % Filename ist in den UserData der TSTOOL-Figure gespeichert
handles = get(fighandle, 'UserData');	% structure with handles to the different menus/controls of the figure
answer=get(handles.loadhandle,'UserData');
TSTOOLfilter=answer{1};
cwd=answer{2};
datafiles=get(handles.lboxhandle,'UserData');
ploteverytime = 0;
if strcmp(get(handles.plotmodehandle(1), 'Checked'), 'on')
  plotmode = 'small';
  ploteverytime = 1;
end

if strcmp(get(handles.plotmodehandle(2), 'Checked'), 'on')
  plotmode = 'large';
  ploteverytime = 1;
end

newfilename = '';
sig = signal(filename);
sig2 = signal(secfile);

if ~isempty(sig)
  if ~isempty(sig2)
    newsig = feval(func, sig, sig2, varargin{:});
    
    if ~isempty(newsig)
      newfilename = newname(filename, [func '_'], '', TSTOOLdatapath);
      newsig=write(newsig, newfilename);
      if ~isempty(newfilename)
	[status,datafiles]=siglevel(newsig,datafiles);
	datafiles=sortdatafiles(datafiles);
	filllistbox(datafiles,handles,name(newsig));
	if ploteverytime
%	  if strcmp(plotmode, 'small')
%	    view(newsig, 6);
%	    %set(gca, 'Visible', 'off');
%	    zoom off
%	    rotate3d off
%	  else
	    tsplot(newfilename, plotmode)
%	  end
	end
      end
    else
      return;
    end
    
  end
end


