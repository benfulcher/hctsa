function mgui(action)
% MGUI   GUI for data analysis programmes.
%    MGUI starts a GUI and supplies Matlab programmes
%    for their easy application to data which are in 
%    the Matlab workspace. 
%
%    The presented programmes are stored in the file
%    MGUI.RC where own programmes can be added. Just
%    include a line with the name of the programme and
%    the minimal and maximal number of arguments, divided
%    by a blank space or tabulator as a separator, e.g.
%
%       plot 1 4
%
%    If the embedded programmes provide an output, then 
%    this output will be stored in the variable ANS in
%    the Matlab workspace.
%
%    See also MGUICLEAN.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2002-2008
% Christian Hoennicke, Norbert Marwan, Andre Sitz, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
% 
%     * Redistributions of source code must retain the above
%       copyright notice, this list of conditions and the following
%       disclaimer.
%     * Redistributions in binary form must reproduce the above
%       copyright notice, this list of conditions and the following
%       disclaimer in the documentation and/or other materials provided
%       with the distribution.
%     * All advertising materials mentioning features or use of this 
%       software must display the following acknowledgement:
%       This product includes software developed by the University of
%       Potsdam, Germany, the Potsdam Institute for Climate Impact
%       Research (PIK), and its contributors.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
% CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
% BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS
%
% $Date: 2010/01/06 08:59:52 $
% $Revision: 1.14 $
%
% $Log: mgui.m,v $
% Revision 1.14  2010/01/06 08:59:52  marwan
% change from GPL to BSD license
%
% Revision 1.13  2009/03/24 08:32:57  marwan
% copyright address changed
%
% Revision 1.12  2006/03/29 13:07:55  marwan
% problems regarding OPRPs and embedding resolved
%
% Revision 1.11  2006/03/16 14:55:57  marwan
% code flattened
%
% Revision 1.10  2006/02/14 11:45:49  marwan
% *** empty log message ***
%
% Revision 1.9  2004/11/10 07:05:21  marwan
% initial import

warning off

if nargin==0
  action='init';
end

    font.FontName=get(0,'factoryTextFontName');
    font.FontUnits=get(0,'factoryTextFontUnits');
    font.FontSize=get(0,'factoryTextFontSize');
    font.FontWeight=get(0,'factoryTextFontWeight');
    font.FontAngle=get(0,'factoryTextFontAngle');

    props.root.defaultUIControlBackgroundColor=get(0,'factoryUIControlBackgroundColor');
    props.root.defaultUIControlForegroundColor=get(0,'factoryUIControlForegroundColor');
    props.root.defaultTextColor=get(0,'factoryUIControlForegroundColor');
    props.window.Color=[0.8008 0.7500 0.6875];
    props.window.Units='Char';
    props.window.Colormap=[[0 0 0];[props.window.Color]];
    props.newwindow.Color=[.8 .8 .8];
    props.newwindow.Units='Pixels';
    props.frame=font;
    props.frame.Style='frame';    
    props.frame.BackgroundColor=props.window.Color;
    props.frame.ForegroundColor=[0 0 0];
    props.frame.Units='Char';
    props.listbox=font;
    props.listbox.FontName='lucida';
    props.listbox.FontSize=10;
    props.listbox.Style='listbox';    
    props.listbox.BackgroundColor=[0.9297 0.8711 0.7969];
    props.listbox.ForegroundColor=[0 0 0];
    props.listbox.Units='Char';
    props.button=font;
    props.button.Style='pushbutton';    
    props.button.BackgroundColor=[0.8008 0.7500 0.6875];
    props.button.ForegroundColor=[0 0 0];
    props.button.Units='Char';
    props.checkbox=font;
    props.checkbox.Style='checkbox';
    props.checkbox.HorizontalAlignment='left';
    props.checkbox.BackgroundColor=props.window.Color;
    props.checkbox.ForegroundColor=[0 0 0];
    props.checkbox.Units='Char';
    props.text=font;
    props.text.Style='text';
    props.text.HorizontalAlignment='left';
    props.text.BackgroundColor=props.window.Color;
    props.text.ForegroundColor=[0 0 0];
    props.text.Units='Char';
    props.msgbox=font;
    props.msgbox.FontSize=get(0,'factoryUIControlFontSize');
    props.msgbox.BackgroundColor=props.window.Color;
    props.msgbox.ForegroundColor=[0 0 0];
    props.msgboxwin.Color=props.msgbox.BackgroundColor;
    v=version;
    if str2double(v(findstr(v,'(R')+2:findstr(v,')')-1))<= 12;
       props.msgboxwin.Colormap=[[0 0 0];[props.msgbox.BackgroundColor]];
    end
    props.logo.Visible='off';
    props.logo.Units='Char';
    props.logo.Color=props.window.Color;

    set(0,'Units','Character'), sc=get(0,'ScreenSize'); 
    set(0,'Units','Pixels'), sc2=get(0,'ScreenSize');
    po=get(0,'defaultFigurePosition');
    props.window.Position=[4 sc(4)-32 125.1667 29.0714]; % Mainframe
    props.newwindow.Position=[sc2(3)-1.05*po(3) po(2) po(3) po(4)]; % Ausgabefenster
    props.frame.Position=[28.0000  5.5 94.8333 21.67]; % frame um die parameterliste
    props.listbox_function_select.Position=[  3.1667  5.5 23.5000 22.1000]; % Funktionselect [.1 .8 .8 .1]
    props.button_apply.Position=[ 99.8333  1.8571 16.8333  2.1429]; % Knopf 1 :: Ausführen [.52 .05 .38 .10]
    props.button_parm_le.Position=[ 28.9 24.3  3.5000  1.5000]; % parm_le :: Parameter verringern [.1 .1 .3 .1]
    props.button_parm_ge.Position=[117.9 24.3  3.5000  1.5000]; % parm_ge :: Parameter erhöhen [.1 .2 .3 .1]
    props.listbox_parm.Position=[ 34.  7.9000 18.5000 16.5000]; % PopupBox (Listbox)
    props.button_ok.Position=[ 34.  6.4 18.5000  1.5000]; % Knopf schliesen
    props.button_parm.Position=[ 34. 24.3 18.5000  1.5000]; % ParameterButton
    props.button_help.Position=[ 76.5000  1.8571 10.1667  2.1429]; % Knopf 5 :: Help
    props.button_close.Position=[ 88.1667  1.8571 10.1667  2.1429]; % Knopf 6 :: Close
    props.text_disclaimer.Position=[  11.5  0.7 23.5000  3.6429]; % Disclaimer
    props.logo.Position=[  3.1667  0.9 6.8 3.5]; % Disclaimer
    props.checkbox_CheckNewWindow.Position=[  36.1667  1.3571 35.1667  1.5000]; % CheckNewWindow
    props.checkbox_ForceSameLength.Position=[ 36.1667  3.0000 35.1667  1.5000]; % ForceSameLength
    props.text_funktion_helptext.Position=[ 29 26.2857 93.5000  1.2857]; % funktion_helptext

set(0,'defaultUIControlBackgroundColor',props.msgbox.BackgroundColor,...
      'defaultUIControlForegroundColor',props.msgbox.ForegroundColor,...
      'defaultTextColor',props.text.ForegroundColor)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splash the BSD

filename='mgui';
which_res=which([filename,'.m']);
bsdrc_path=[strrep(which_res,[filename,'.m'],''), 'private'];
bsdrc_file=[bsdrc_path, filesep, '.bsd.',filename];
if ~exist(bsdrc_path,'dir')
  mkdir(strrep(which_res,[filename,'.m'],''),'private')
end
if ~exist(bsdrc_file,'file')
  fid=fopen(bsdrc_file,'w');
  fprintf(fid,'%s\n','If you delete this file, the BSD License will');
  fprintf(fid,'%s','splash up at the next time the programme starts.');
  fclose(fid);

  if exist('bsd')
    txt=bsd;
  else
    txt={'The BSD License was not found.'};
  end
  h=figure('NumberTitle','off',...,
         'ButtonDownFcn','close',...
         'Name','BSD License');
  ha=get(h,'Position');
  h=uicontrol('Style','Listbox',...
            'ButtonDownFcn','close',...
            'CallBack','close',...
            'Position',[0 0 ha(3) ha(4)],...
            'FontName','Courier',...
            'BackgroundColor',[.8 .8 .8],...
            'String',txt);
  waitfor(h)
end


hMF = findobj('Tag','MainFrame');
if ~(isempty(hMF)), props = get(hMF,'UserData'); end
if ~isempty(findobj('Tag','helpdlgbox')), delete(findobj('tag','helpdlgbox')), end
hfs    = findobj('Tag','function_select');
select = get(hfs,'Value');
funktion = get(hfs,'UserData');
if ~isempty(funktion)
  FktName = char(funktion(select).name);
  if strcmp(FktName,'------------');
    select=select+1;
    set(hfs,'Value',select)
    mgui action
  end
end
set(0,'ShowHiddenHandles','off')

% Und "Aktion!" - hier geht es los
switch(action)

  case 'init'

  %%% Hier wird der Programmablauf vorbereitet,
  %%% * die Graphikkoordinaten festgelegt
  %%% * die Funktionslisten aus der Datei generiert
  %%% * das Grundfenster gezeichnet

   % Lösche bestehende Fenster gleichen Tags
    if ~isempty(findobj('Tag','MainFrame'))
      delete(findobj('Tag','MainFrame'))
    end

    % Die Funktionen Load und Save am Anfang einfuegen
    i = 1;
    funktion(i).name = cellstr('load');
    funktion(i).min = 0;
    funktion(i).max = 0;
    funktion(i).akt = 0; 
    funktion(i).klicklist{4}=[];
    funktion(i).klickstring{4}=[];
    i=i+1;

    funktion(i).name = cellstr('save');
    funktion(i).min = 1;
    funktion(i).max = 1;
    funktion(i).akt = 1; 
    funktion(i).klicklist{4}=[];
    funktion(i).klickstring{4}=[];
    i=i+1;

    funktion(i).name = cellstr('add programme');
    funktion(i).min = 0;
    funktion(i).max = 0;
    funktion(i).akt = 0; 
    funktion(i).klicklist{4}=[];
    funktion(i).klickstring{4}=[];
    i=i+1;

    funktion(i).name = cellstr('------------');
    funktion(i).min = 0;
    funktion(i).max = 0;
    funktion(i).akt = 0; 
    funktion(i).klicklist{4}=[];
    funktion(i).klickstring{4}=[];

    %% Einlesen der Funktionsliste
    % öffnen der Datei mit den Funktionensnamen und -parametern
    fid = fopen('mgui.rc','r');
    if (fid ~= -1) % fopen liefert im Fehlerfall "-1"
      line = fgetl(fid); 
      while (line ~=-1) %EOF
        if ~isempty(line) & line(1)~='%' %einlesen erster Parameter
          i = i + 1;
          [name,line] = strtok(line);
          funktion(i).name = cellstr(name);
          [min,line] = strtok(line);
          if ~isempty(min)
            funktion(i).min = str2double(min);
            max = strtok(line);
            if ~isempty(max) 
              % Hack : mehr als 4 Parameter werden nicht zugelassen
              funktion(i).max = str2double(max);
              if (funktion(i).max > 4)
                funktion(i).max = 4;
              end
              funktion(i).akt = str2double(min);
              funktion(i).klicklist{4}=[];
              funktion(i).klickstring{4}=[];
              % Fehlerfall min > max
              if funktion(i).min > funktion(i).max
                h=errordlg([' Error in data structure at: ' char(funktion(i).name) ' -> (min>max)'],'Initialization');
                set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
                waitfor(h)
              end
         % Fehler: "Zu wenige Parameter" einsammeln
            else
              h=errordlg([' Error in data structure at: ' char(funktion(i).name) ' -> (to less parameter)'],'Initialization');
              set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
              waitfor(h)
              funktion(i).min = 0;
              funktion(i).max = 0;
              funktion(i).akt = 0;
              funktion(i).klicklist{4}=[];
              funktion(i).klickstring{4}=[];
            end
          else
            h=errordlg([' Error in data structure at: ' char(funktion(i).name) ' -> (to less parameter)'],'Initialization'); 
            set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
            waitfor(h)
            funktion(i).min = 0;
            funktion(i).max = 0;
            funktion(i).akt = 0;
            funktion(i).klicklist{4}=[];
            funktion(i).klickstring{4}=[];
          end
        end
        line = fgetl(fid); % Zeilenweises lesen
      end
      fclose(fid); % schliesen der Datei

    % Wenn fopen() fehlgeschlagen ist
    else
      error('File mgui.rc not found.');
%       funktion.name = cellstr('Empty');
%       funktion.min = 0;
%       funktion.max = 0;
%       funktion.akt = 0;
%       funktion.klicklist{4}=[];
%       funktion(i).klickstring{4}=[];
    end

    % Schmeiße alle unbekannten Funktionen raus...
    j = size(funktion,2); i=5;
    while i<=j
      if isempty(double(which(char(funktion(i).name))))
        funktion(i) = [];
        i = i-1;
        j = j-1;
      end
        i=i+1;
    end

    %% Zeichnen der Obefläche
    % Zeichnen des Grundfensters (Hack: Koordinaten sind fix)
    figure(props.window,'Tag','MainFrame',...
           'NumberTitle','off',...
           'Name','Controlpanel',...
           'MenuBar','None',...
           'DeleteFcn','mgui winclose',...
           'UserData',props,...
           'Resize','Off');

    uicontrol(props.frame);
    
    if length(funktion)>4, select=5; else select=3; end

    uicontrol(props.listbox, props.listbox_function_select,...
              'Tag','function_select',... 
              'String', [funktion.name]',...
              'UserData', funktion,...
              'Value',select,...
              'Callback', 'mgui mod2', 'ButtonDownFcn','helpwin mgui');

    % Knopf 3 :: Parameter verringern
b=[fliplr(triu(ones(10,10))),triu(ones(10,10))];
b=rot90(b(:,1:2:end));
for i=1:3,
c=(b==1)*props.button.BackgroundColor(i)+(b==0)*props.button.ForegroundColor(i);
icon(:,:,i)=c(:,1:2:end);
end
    uicontrol(props.button, props.button_parm_le,...
              'Tag','parm_le',...
              'CData',icon,...
              'Visible', 'off',...
              'Callback','mgui ChangeNumberOfParameters',... 
              'ToolTip','Decreases the number of parameters.');

    for i =1:4,
        uicontrol(props.button, 'Tag',['button_parm' num2str(i)],...
                  'Position', [(props.button_parm.Position(1) + ((i - 1) * 21.2)), props.button_parm.Position(2), props.button_parm.Position(3), props.button_parm.Position(4)],...
                  'String', ['Parm' num2str(i)],...
                  'Visible', 'off',...
                  'Callback', 'mgui Vis');
                  
        uicontrol(props.listbox, 'Tag',['listbox_parm' num2str(i)],...
                  'Value',[],...
                  'Max',2,...  %Hack: scheinbar beliebig Zahl(>1) möglich...
                  'ButtonDownFcn','mgui selectcolumn',...
                  'Visible', 'off',...
                  'Position', [(props.listbox_parm.Position(1) + ((i - 1) * 21.2)), props.listbox_parm.Position(2), props.listbox_parm.Position(3), props.listbox_parm.Position(4)],...
                  'Callback','mgui klickreihenfolge'); 
                  
        uicontrol(props.button, 'Tag',['button_close_parm' num2str(i)],...
                  'String', 'OK',...
                  'Visible', 'off',...
                  'Position', [(props.button_ok.Position(1) + ((i - 1) * 21.2)), props.button_ok.Position(2), props.button_ok.Position(3), props.button_ok.Position(4)],...
                  'Callback', 'mgui UnVis',...
                  'UserData',[]);
    end


    % Knopf 4:: Parameter erhöhen
 b=[fliplr(triu(ones(10,10))),triu(ones(10,10))];
 b=rot90(b(:,1:2:end),-1);
 for i=1:3,
 c=(b==1)*props.button.BackgroundColor(i)+(b==0)*props.button.ForegroundColor(i);
 icon(:,:,i)=c(:,1:2:end);
 end

    uicontrol(props.button, props.button_parm_ge,...
              'Tag','parm_ge',...
              'CData',icon,...
              'Visible', 'on',...
              'Callback','mgui ChangeNumberOfParameters',... 
              'ToolTip','Increases the number of parameters.');

    uicontrol(props.checkbox, props.checkbox_ForceSameLength, ...
              'Tag','ForceSameLength',...
              'Value',0,...
              'ToolTip','Fills parameter vectors with NaN in order to get the same vector lengths.',...
              'String','Force equal data length.');

    uicontrol(props.checkbox, props.checkbox_CheckNewWindow, ...
              'Tag','CheckNewWindow',...
              'Value',0,...
              'ToolTip','Every apply click will create a new figure.',...
              'String','Create new figure.');

    % Knopf 5 :: Help
    uicontrol(props.button, props.button_help,...
              'Tag','help',...   
              'String','Help',...
              'Callback', 'mgui help',...
              'ToolTip','Helpwindow.', 'ButtonDownFcn','helpwin mgui');

    % Knopf 6 :: Close
    uicontrol(props.button, props.button_close, ...
              'Tag','close',...
              'String','Close',...
              'Callback', 'mgui close',...
              'ToolTip','Closing the GUI.', 'ButtonDownFcn','helpwin mgui');

    % Knopf 1 :: Ausführen
    uicontrol(props.button, props.button_apply, 'Tag','apply',...  
              'String','Apply',...
              'Enable','off',...
              'Callback', 'mgui apply',...
              'ToolTip','Apply the choosen function.', 'ButtonDownFcn','helpwin mgui');

    h=axes;
    logo=load('logo');
    imagesc([logo.logo fliplr(logo.logo)], 'ButtonDownFcn','helpwin mgui')
    set(h,props.logo)

    h=uicontrol(props.text, props.text_disclaimer, ...
              'HandleVis','off',...
              'Tag','disclaimer',...
               'ButtonDownFcn','helpwin mgui');
    h2=textwrap(h,{[char(169),' AGNLD'],'University of Potsdam','2002-2006'});
    set(h,'String',h2)
    
 
     mgui mod2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'mod2'

    %%% erlegt den ganzen init-kram, wenn eine Funktion ausgewählt wurde

    % Anzahl der Parameter und Auswahlboxen beim Start anpassen und zeichnen
    % lesen der Anzahl der minimalen/maximalen Parameter
    FN_Min  = funktion(select).min;
    FN_Max  = funktion(select).max;
    FN_Akt  = funktion(select).akt;
    FN_Name = funktion(select).name;

%%% kurzhilfe einblenden
    if strcmp(FN_Name, 'add programme');
      fnhelptext='ADD PROGRAMME   Adds a programme into the MGUI and its related RC file.';
      ind=length(fnhelptext);
    else
      fnhelptext=help(char(FN_Name));
      ind=findstr(char(10),fnhelptext);
    end

    hfh=findobj('Tag','funktion_helptext');
    if isempty(hfh)
      hfh=uicontrol(props.text, props.text_funktion_helptext,...
                'Tag','funktion_helptext',...
                'String',fnhelptext(1:ind(1)),...
                'ButtonDownFcn','helpwin mgui');
    else
      set(hfh,'string',fnhelptext(1:ind(1)))
    end
    h=get(hfh,'Extent');h2=get(hfh,'Position');
    set(hfh,'Position',[h2(1:2) h(3) h2(4)])
    drawnow

    %% Knöppe in Ausgangsstellung bringen
    % Wenn, mehrere Parameter erlaubt sind, die Knöpfe aktivieren
    % Es sind mehr Parameter möglich

    set(findobj('Tag','parm_le'),'Visible','on');
    set(findobj('Tag','parm_ge'),'Visible','on');

    if (FN_Akt == FN_Max)
      set(findobj('Tag','parm_ge'),'Visible','off');
    end

    if (FN_Akt == FN_Min)
      set(findobj('Tag','parm_le'),'Visible','off');
    end

    check_workspace;
    check_apply

    %% aufdecken 
    for i = 1 : FN_Akt
      hparm = findobj('Tag',['button_parm' num2str(i)]);
      hlbox = findobj('Tag',['listbox_parm' num2str(i)]);
      hbcp  = findobj('Tag',['button_close_parm' num2str(i)]);
      set(hparm,'Visible','on');
      set(hlbox,'Visible','on');
      set(hbcp,'Visible','on');
    end

    %% Verdeckt
    for i = (FN_Akt +1) : 4 % HACK: es sind nicht mehr als 4 paramter vorgesehen
      hparm = findobj('Tag',['button_parm' num2str(i)]);
      hlbox = findobj('Tag',['listbox_parm' num2str(i)]);
      hbcp  = findobj('Tag',['button_close_parm' num2str(i)]);
      set(hparm,'Visible','off');
      set(hlbox,'Visible','off');
      set(hbcp,'Visible','off');
    end

    if strcmp(FN_Name, 'load');
       mgui load
    end

    if strcmp(FN_Name, 'add programme');
       mgui add_programme
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
  case 'apply' 
    %%% Alle Parameter werden ausgelesen und hintereinander zusammengestellt und ausgeführt
             
    %% Init
    FN_Akt = funktion(select).akt;
    FktName = char(funktion(select).name);

    %% 'save' und 'load' rausfischen -> selbst geschrieben
    if strcmp(FktName,'save')
      mgui save;
      return;
    end 

    if strcmp(FktName,'load')
      mgui load;
      return;
    end 

    if strcmp(FktName, 'add programme');
      mgui add_programme;
      return;
    end

    %% Auszuführenden String zusammenbasteln
    Argumente = '(';
    for j=1 : FN_Akt
      Arg = ['Arg',num2str(j)];
%    %% Init
      var_list = get(findobj('Tag',['listbox_parm' num2str(j)]),'UserData');
      temp=funktion(select).klicklist;
      if ~isempty(temp), klick_list=temp{j}; else klick_list=[]; end
      if isempty(klick_list), klick_list=get(findobj('Tag',['button_close_parm', num2str(j)]),'UserData'); end

      Matrix = []; % Zu speichernder Inhalt

      %% Parameterinhalt zusammenstellen
      for i = 1:size(klick_list,2); 
        % Einlesen des Variablennamens und zuweisen seines Inhaltes
        VarName = var_list(klick_list(1,i)).name;

        try
          VarValue = evalin('base',VarName);
        catch
          h=warndlg(['Variable ' VarName ' could not be found and will be ignored henceforth.'],'Scan workspace'); 
          set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
          set(h,'Tag','helpdlgbox','HandleVisibility','On');
          continue;
        end

        % Zeilenvektor in Spaltenvektor transformieren 
        if (size(VarValue,1)) == 1
          VarValue = VarValue';
        end

        % Feststellen der gewünschten Reihen
        SelectRows = klick_list(klick_list(:,i) ~= 0,i);
        SelectRows(1,:) = []; % Lösche Kopf (Enthält Verweis auf den Variablennamen)

        % Wenn Spalten gewählt wurden, dann den Varibleninhalt darauf begrenzen
        if ~isempty(SelectRows)
           VarValue = VarValue(:,SelectRows);
        end

        % Unterschiedliche Länge -> Auffüllen mit "NaN"
        if (size(Matrix,1) ~= size(VarValue,1)) & (size(Matrix,1) ~= 0)
          % Matrix zu lang
          if (size(Matrix,1) > size(VarValue,1))  
             VarValue((size(VarValue,1)+1):size(Matrix,1),:) = NaN;
          % Spaltenvektor zu lang
          else
            Matrix((size(Matrix,1)+1):size(VarValue,1),:) = NaN;
          end
        end

        Matrix = [Matrix VarValue]; 
      end


    % zusammenfügen aller gewählten Variablen

        % Wertzuweisung eines Paramters
        eval([Arg '= Matrix;']);

        Argumente = strcat(Argumente, Arg);

%        if j < FN_Akt
        Argumente = strcat(Argumente,',');    
%        end    
    end
    Argumente(end) = [];
    Argumente = strcat(Argumente,')');

    %% Alle Parameter erhalten die gleiche Länge (wenn gewünscht)
    if get(findobj('Tag','ForceSameLength'),'Value')

      % Init
      MaxLength = 0;

      % Suche des längsten Parameters
      for i = 1:FN_Akt
        Matrix = eval(['Arg' num2str(i)]);
        if (size(Matrix, 1) > MaxLength)
          MaxLength = size(Matrix, 1);
        end
      end

      % Alle Parameter auf die gleiche Länge trimmen
      for i = 1:FN_Akt
        Matrix = eval(['Arg' num2str(i)]); % laden
        Matrix(size(Matrix,1):MaxLength,:)= NaN;
        eval([['Arg' num2str(i)] '= Matrix;']); % speichern
      end

    end

    %% Ausgabe der Funktionswerte
    % Erstellt auf Wunsch ein neues Ausgabefenster
    hnpw = findobj('Tag','CheckNewWindow');
    bNewWindow = get(hnpw,'Value');
    if bNewWindow
        figure(props.newwindow,'Tag','PrintWindow');
    end

    %% Ausführen der ganzen Mühe (mit Errorhandling)

    % Schreibschutz aktivieren 
    set(hMF, 'HandleVisibility','Off');
    % Ausführen 

    try
      set(0,props.root);
      eval([FktName Argumente]);
%       if ~isempty(ans)
%         assignin('base','ans',ans)
%       end
    catch
      set(0,'defaultUIControlBackgroundColor',props.msgbox.BackgroundColor,...
      'defaultUIControlForegroundColor',props.msgbox.ForegroundColor,...
      'defaultTextColor',props.text.ForegroundColor)
      h=errordlg(lasterr,['Call ',upper(FktName)]);
      set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
      waitfor(h);
    end

    % Schreibschutz deaktivieren
    set(hMF, 'HandleVisibility','on');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'ChangeNumberOfParameters'

    %%% verändern der aktuellen Parameteranzahl
    %%% parm_le => verringern; parm_ge => vergroessern

    %% Init
    hcall  = gcbo;
    caller = get(hcall,'Tag');

    FN_Akt = funktion(select).akt;
    FN_Min = funktion(select).min;
    FN_Max = funktion(select).max;
    
    %% Anzahl sichtbarer Parameter verringern
    if strcmp('parm_le',caller)
      if FN_Akt > FN_Min 
        set(findobj('Tag',['button_parm' num2str(FN_Akt)]),'Visible','off');
        set(findobj('Tag',['button_close_parm' num2str(FN_Akt)]),'Visible','off');
        set(findobj('Tag',['listbox_parm' num2str(FN_Akt)]),'Visible','off');

        funktion(select).akt = FN_Akt - 1;
        set(hfs,'UserData',funktion);
        if funktion(select).akt <= FN_Min
          set(hcall,'Visible','off');
        end
        set(findobj('Tag','parm_ge'),'Visible','on');
      end

      %if ~isempty(get(findobj('Tag',['button_close_parm', num2str(FN_Akt-1)]),'UserData'))
      %  set(findobj('Tag','apply'),'Enable','on');
      %end
    %% Anzahl sichtbarer Parameter erhöhen (parm_ge)
    else
      if FN_Akt < FN_Max 
        funktion(select).akt = FN_Akt + 1;

        set(findobj('Tag',['button_parm' num2str(FN_Akt+1)]),'Visible','on');
        set(findobj('Tag',['button_close_parm' num2str(FN_Akt+1)]),'Visible','on');
        set(findobj('Tag',['listbox_parm' num2str(FN_Akt+1)]),'Visible','on');

        set(hfs,'UserData',funktion);
        if funktion(select).akt >= FN_Max 
          set(hcall,'Visible','off');
        end
        set(findobj('Tag','parm_le'),'Visible','on');
        % Den Apply-Knopf deaktivieren
       % if isempty(get(findobj('Tag',['button_close_parm', num2str(FN_Akt+1)]),'UserData'))
       %   set(findobj('Tag','apply'),'Enable','off');
       % end
      end
    end

    %% Alle "Klick_listen" durchsuchen, ob "Apply" ausgeführt werden darf
    FN_Akt   = funktion(select).akt;

    trigger = 0;
    % Alle "Klick_list"en durchsuchen
    for i = 1 : FN_Akt
      hbcp = findobj('Tag',['button_close_parm' num2str(i)]);
      temp=funktion(select).klicklist;
      if ~isempty(temp), klick_list=temp{i}; else klick_list=[]; end
      if isempty(klick_list), klick_list=get(hbcp,'UserData'); end
      if isempty(klick_list)
        trigger = 1;
      end
    end
    % Wenn mindestens eine leere dabei, dann die Ausführung verhindern
    if trigger
      set(findobj('Tag','apply'),'Enable','off');
    else
      set(findobj('Tag','apply'),'Enable','on');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
  case 'Vis'
    %%% Macht die geforderte Listbox + Closebutton sichtbar

    % Aufrufender = akt. Parameter wird rot gefärbt
    hcall = gcbo;

    %% Bestimmen des zugehörigen Auswahlfensters und Closebuttons
    caller = get(hcall,'Tag');
    caller_lbox = strrep(caller,'button','listbox');
    caller_bc   = strrep(caller,'button','button_close');

    lokal_ws=check_workspace;

%%%%%%%%%%%%%%
% check for erased variables
nparm = str2double(strrep(caller,'button_parm',''));                             
tmplist = funktion(select).klicklist{nparm};                                  
tmpstr = funktion(select).klickstring{nparm};                                 
if ~isempty(tmpstr)                                                           
ParmBereiche =[0 find(tmpstr == ',') size(tmpstr,2)+1];                       
                                                                             
% jeden Parameter einzeln checken                                             
for i = 1:size(ParmBereiche,2)-1                                              
  varname = tmpstr(ParmBereiche(i)+1:ParmBereiche(i+1)-1);                    
  trigger = 0;                                                                
  for j = 1:size(lokal_ws,2)                                                  
    if strcmp(lokal_ws(j).name,varname)                                       
      trigger =1;                                                             
      tmplist(1,i) = j;                                                       
      break;                                                                  
    end                                                                       
  end                                                                         
  if trigger == 0                                                             
    tmplist(i) = NaN;                                                         
    tmpstr(ParmBereiche(i)+1:ParmBereiche(i+1)-1) = ' ';                      
  end                                                                         
end                                                                           
tmplist(isnan(tmplist(1,:))) = [];                                      
end                                                                           
%% hier noch die richtigen pfade eintragen                                    
%set(findobj('tag',caller_bc),'UserData',tmplist);                            
%set(hparm,'String',tmpstr);                                                  
funktion(select).klicklist{nparm}= tmplist;                                   
funktion(select).klickstring{nparm}=tmpstr;                                   
set(hfs,'UserData',funktion);                                                 
%%%%%%%%%%%%%%


    %% Aktivieren der angeforderten Listbox + Closebutton
    % Erst Listbox
    hlbox=findobj('Tag',caller_lbox);
    set(hlbox,'Visible','on');
%    set(hlbox,'String',[Ausgabe]);
%    set(hlbox,'UserData',lokal_ws);
    if isempty(get(hcall,'tooltip'))
      set(hlbox,'Value',[]); % beim oeffnen ist nichts ausgewaehlt
    end

    % Dann Closebutton
    hbcp=findobj('Tag',caller_bc);
    set(hbcp,'Visible','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'UnVis'
    %%% macht die Listbox + Closebutton unsichtbar

    %% bestimmen des zugehörigen Auswahlfensters und Closebuttons
    hcall = gcbo;
    caller = get(hcall,'Tag');
    caller_lbox = strrep(caller,'button_close','listbox');

    %% deaktivieren der angeforderten Elemente
    hlbox=findobj('Tag',caller_lbox);
    set(hlbox,'Visible','off');
    set(hcall,'Visible','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'help'
    %%% Es werden die die Standart-Matlab-Hilfetexte angezeigt

    %% Entsprechenden Hilfetext auswählen und anzeigen
    if ~isempty(select)
      % 'Load' und 'Save' haben eigenen Hilfstext
      if strcmp(funktion(select).name,'load') | strcmp(funktion(select).name,'save') | strcmp(funktion(select).name,'add programme')
        if strcmp(funktion(select).name,'load')
          h=helpdlg(['Load workspace variables from an extern data file.',char(10),'Additionally text lines (e.g. header lines) will be skipped.'],'Help');
        elseif strcmp(funktion(select).name,'save') % (funktion(select).name == 'SAVE')
          h=helpdlg(['Save workspace variables into an extern data file.',char(10),'The file format will be plain ASCII.'],'Help');
        else
          h=helpdlg(['Adds a new programme into the function list of this MGUI.',char(10),...
                     'The minimal and maximal number of input parameters for this',char(10),...
                     'programme has to be specified. Finally, this changing will', char(10),...
                     'be done in the mgui.rc file.'],'Help');
        end
        set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
         set(h,'Tag','helpdlgbox','HandleVisibility','On');
   else
      % Hilfstext
        eval(['helpwin ', char(funktion(select).name)]);
      end
    else % wenn keine Funktion ausgewählt wurde
      helpwin mgui;
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'close'
    % beendet das Programm
    close(hMF)
    warning on
    % clear all
    tx{1}='Thank you for your blind trust.';
    tx{2}='Do you really like this?';
    tx{3}='Thanks, that you believe in us.';
    tx{4}='We hope we could satisfy you.';
    tx{5}='Thank you for working with us.';
    tx{6}='Now you should go home.';
    tx{7}='Did you already find us in the WorldWideWeb?';
    tx{8}='Your cat might be hungry.';
    tx{9}='Maybe your Matlab license will be needed by others?';
    tx{10}='Can you imagine, how long we have worked for that?';
    tx{11}='Thank you for your responsiveness.';
    tx{12}='Now ready for publishing?';
    disp(tx{round((length(tx)-1)*rand)+1})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'winclose'
    % koennte noetig werden, wenn noch irgendwelche settings 
    % zurueckgestellt werden muessen
    set(0,props.root);
    warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'klickreihenfolge'

    %%% Die Klickreihenfolge wird im parm_button_close gespeichert

    %% Bestimmung des Aufrufers
    hlbox = gcbo;

    % die passenden Handel initialisieren
    tmp = strrep(get(hlbox,'Tag'),'listbox','button_close');
    hbcp   = findobj('Tag',tmp);

    tmp = strrep(get(hlbox,'Tag'),'listbox','button');
    hparm = findobj('Tag',tmp);

    %% Init
    nparm = str2double(strrep(get(gcbo,'Tag'),'listbox_parm',''));
    temp=funktion(select).klicklist;
    if ~isempty(temp), 
       klick_list=temp{nparm}; 
    else
       klick_list = [];
    end
    if isempty(klick_list), klick_list = get(hbcp,'UserData'); end
    klick_akt  = get(hlbox,'Value');
    if isempty(klick_akt), return; end; % Keine Variable gewählt

    var_list   = get(hlbox,'UserData');
    if isempty(var_list(1).name), return; end; % Keine Variable vorhanden

    %% Akualisieren der Klick_list (und anzeigen)
    % Wenn ausversehen mehr als ein Wert gewählt wurde
    if (size(klick_akt,2) > 1)
      h=errordlg('Choose only one variable.','Parameter selection');
      set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
      waitfor(h)
      set(hlbox,'Value',klick_list(1,:));
      return;
    else

    % Aktualisieren
    if isempty(klick_list)
      posi=[];
    else
      posi = find(klick_list(1,:) == klick_akt);
    end

    if(isempty(posi))
      % Einfügen
      klick_list(1,(end +1)) = klick_akt;
    else
      % Löschen
      klick_list(:,posi) = [];
    end

    % Das der Nutzer sein klicken sieht (Listbox)
    set(hlbox,'Value',klick_list(1,:));

    klick_string=make_klick_string(klick_list,var_list);

       set(hparm,'String',klick_string);
       set(hparm,'ToolTip',klick_string);

      % Speichern der ganzen Mühe
      set(hbcp,'UserData',klick_list);
      funktion(select).klicklist{nparm}=klick_list;
      funktion(select).klickstring{nparm}=klick_string;
      set(hfs,'UserData',funktion)
    end
     
    check_apply

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case('selectcolumn')
    
    %%% Hier erfolgt die Auswahl der Spalten

    %% Bestimmung des Aufrufers
    hcall = gcbo;
    caller= get(hcall,'Tag');

    % die passenden Handel initialisieren
    tmp   = strrep(caller,'listbox','button_close');
    hbcp   = findobj('Tag',tmp);

    hlbox = hcall;

    tmp = strrep(caller,'listbox','button');
    hparm = findobj('Tag',tmp);

    %% Init
    nparm = str2double(strrep(get(hlbox,'Tag'),'listbox_parm',''));
    temp=funktion(select).klicklist;
    if ~isempty(temp), klick_list=temp{nparm}; else klick_list = []; end
    if isempty(klick_list), klick_list = get(hbcp,'UserData'); end
    if isempty(klick_list), return; end; % Keine Variable ausgewählt
    var_list   = get(hlbox,'UserData');
    if isempty(var_list(1).name), return; end; % Keine Variable vorhanden

    klick_string = [];
    VarName = var_list(klick_list(1,end)).name;
    VarValue = evalin('base',VarName);
    [Zeilen Spalten] = size(VarValue); % Begrenzungen ermitteln  

    %% Einlesen der gewünschten Spalten   
    % Zeilenvektor -> wie Spaltenvektor behandeln
    if Zeilen == 1
      [Spalten Zeilen] = size(VarValue); % Vertauschen 
    end

    % Wenn vorhanden, dann soll die alte Auswahl angezeigt werden
    oldcolumns = klick_list(2:end,end);
    oldcolumns = oldcolumns(oldcolumns(:) ~= 0); % Nur gewähltes zeigen
    oldcolumns = cellstr(num2str(oldcolumns')); % Umwandeln zu Textausgabe
    if isempty(oldcolumns)
      oldcolumns = cellstr(' ');
    end

    % Eingabe der Spalten
    set(0,'DefaultUIControlBackgroundColor',props.msgbox.BackgroundColor);
    answer = inputdlg('Use '','' or '':'' or '' '' for selecting the columns.','Multi-column vector',1,oldcolumns);
    set(0,props.root);

    % Wenn etwas eingegeben wurde
    if ~isempty(answer) 
      if strcmp(answer,':') | strcmp(answer,' '), answer='0'; end
      answer3 = str2double(char(answer)); % Umwandel in Zahlen
      answer3(answer3 > Spalten) = []; % Gewähltes <= Max. Spalten

      % Fehler sammeln -> Verarbeitung abbrechen
      if isempty(answer3)
        h=errordlg('Wrong input. Please use Matlab notation.','Multi-column vector');
        set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
        waitfor(h)
        return;
      end
      if size(answer3,1) ~= 1
        h=errordlg('The dimension of the input was multi-dimensional.','Multi-column vector');
        set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
        waitfor(h)
        return;
      end

      answer3 = answer3';
      klick_list(2:Zeilen,end) = 0; % Clean
      klick_list(2:(size(answer3,1)+1),end) = answer3;
    end

    %% Anzeigen der Parameter

    for i = 1:size(klick_list,2)
      klick_string = strcat(klick_string,var_list(klick_list(1,i)).name);
      selected = klick_list(klick_list(:,i) ~= 0,i);

      if size(selected,1) > 1
        klick_string = strcat(klick_string, '('); 

        for j = 2:size(selected,1)
          klick_string = strcat(klick_string,num2str(selected(j)));
          klick_string = strcat(klick_string, ','); 
        end
        klick_string(end) = []; 

        klick_string = strcat(klick_string, ')'); 
      end 
 
     klick_string = strcat(klick_string, ','); 
    end
    klick_string(end) = []; 

    % Speichern der Ausgabe
    set(hparm,'String',klick_string);
    set(hparm,'ToolTip',klick_string);


      % Speichern der Mühe
      set(hbcp,'UserData',klick_list);
      funktion(select).klicklist{nparm}=klick_list;
      funktion(select).klickstring{nparm}=klick_string;
      set(hfs,'UserData',funktion)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case ('load')

    %%% Lädt, so weit es geht, den Inhalt einer frei wählbaren Datei
    %%% unter einem frei wählbaren Namen im "base"

    %% Einlesen der gewünschten Datei
    % öffnen der Datei
    setptr(gcf,'watch')
    if isunix
       [Filename, Pathname]=uigetfile('*');
    else
       [Filename, Pathname]=uigetfile('*.*');
    end
    if Pathname == 0 & Filename == 0, setptr(gcf,'arrow'), return; end; % Wenn 'Cancel'

    fid = fopen([Pathname Filename],'r');
    if (fid ~= -1) % fopen liefert im Fehlerfall "-1"
      Matrix = []; % Variableninhalt
      %i = 0;
      Line = fgetl(fid); 
      while (Line ~=-1) %EOF
        if ~isempty(Line) & Line(1)~='%' %einlesen erster Parameter
          %i = i + 1;
          LineX = [];
          while ~isempty(Line)          
             [VarX,Line] = strtok(Line);
             VarX = str2double(VarX);
             if isempty(VarX) %unbekanntest wird als NaN gefiltert
               VarX = NaN;
             end;
             LineX = [LineX, VarX];
          end
          % Im Bedarfsfall die Matrix vergrössern
          if (size(LineX,2) > size(Matrix,2)) & (size(Matrix,1) ~=0)
            Matrix(:,size(Matrix,2)+1:size(LineX,2)) = NaN;
          end

          % Im Bedarfsfall die Eingelesene Zeile vergrössern
          if size(LineX,2) < size(Matrix,2)
            LineX(1,size(LineX,2)+1:size(Matrix,2)) = NaN;
          end
          % Zeile Anhängen
          Matrix = [Matrix; LineX];
        end
      Line = fgetl(fid); % Zeilenweises lesen
      end

      %% freien Variablennamen finden
      VarName_free = 0; % 0 = false // 1 = true

      while VarName_free == 0
        VarName_free = 1;

        % Variablennamen bestimmen
        set(0,'DefaultUIControlBackgroundColor',props.msgbox.BackgroundColor);
        VarName = inputdlg(['File ',Filename,' successfully read! Rows: ' num2str(size(Matrix,1)) ' Columns: ' num2str(size(Matrix,2)),char(10),'Choose a variable name:'],'New data');
        set(0,props.root);
        if isempty(VarName), break; end; % Wenn 'Cancel'

        % einlesen aller vorhandenen VariablenNamen und ihre Feldgrösse
        lokal_ws = evalin('base','whos');
   
        if length(lokal_ws) > 0
          for i=1:length(lokal_ws) % Alle Variablen Durchtesten
            if strcmp(lokal_ws(i).name,char(VarName))
              VarName_free = 0; % Variablenname existiert schon
              h = errordlg('Variable name already exists.','New data');
              set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
              waitfor(h); % wichtig! -> sonst kommt es zum Konflikt mit dem inputdialog
            end
          end 
        end
      end % while VarName_free
      if isempty(VarName), return; end; % Erster Sprung erfolgt nur aus der Whileschleife

      %% Speichern der Eingelesenen Werte
      assignin('base',char(VarName),Matrix);

      lokal_ws_all = evalin('base','whos');

      if length(lokal_ws_all) > 0
        j = 1;
        for i=1:length(lokal_ws_all)
          if strcmp(lokal_ws_all(i).class,'double') | strcmp(lokal_ws_all(i).class,'single') | strcmp(lokal_ws_all(i).class,'int8') | strcmp(lokal_ws_all(i).class,'uint8') | strcmp(lokal_ws_all(i).class,'int16') | strcmp(lokal_ws_all(i).class,'uint16') | strcmp(lokal_ws_all(i).class,'int32') | strcmp(lokal_ws_all(i).class,'uint32')
            Ausgabe(j,1) = cellstr([lokal_ws_all(i).name, ' [', num2str(lokal_ws_all(i).size(1)),'x',num2str(lokal_ws_all(i).size(2)), ']' ]);
            lokal_ws(j) = lokal_ws_all(i);
            j = j+1;
          end
        end 
      % Der Workspace ist leer 
      else
        h=errordlg('No numeric variables in the workspace.','Scan workspace');
        set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
        set(h,'Tag','helpdlgbox','HandleVisibility','On');
        Ausgabe(1,1) = cellstr('Empty');
        lokal_ws(1).name = [];
      end

      for i=1:4;
        hlbox=findobj('Tag',['listbox_parm',num2str(i)]);
        set(hlbox,'String',Ausgabe);
        set(hlbox,'UserData',lokal_ws);
      end
      
    % Datei konnte nicht geöffnet werden.
    else
      h=errordlg('File not found.','IO process');  
      set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
      waitfor(h)
    end
    fclose(fid); % schliesen der Datei
    set(findobj('Tag','apply'),'Enable','On');
    setptr(gcf,'arrow')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case ('save')

    %%% Speichert die gewählten Variablen und Spalten in einer frei wählbaren Datei
    %%% Der zu speichernde Inhalt liegt "nur" im ersten Parameter

    %% Bestimmen des Speicherortes
    setptr(gcf,'watch')
    [Filename, Pathname]=uiputfile('*.*');
    if Pathname == 0 & Filename == 0, setptr(gcf,'arrow'), return; end; % Wenn 'Cancel'

    %% Init
    temp=funktion(select).klicklist;
    if ~isempty(temp), 
       klick_list=temp{1}; 
    else
       klick_list=[];
    end
    if isempty(klick_list), klick_list=get(findobj('Tag','button_close_parm1'),'UserData'); end
    var_list = get(findobj('Tag','listbox_parm1'),'UserData');
    Matrix = []; % Zu speichernder Inhalt

    if isempty(var_list(1).name), return; end; % Wenn keine Variable vorhanden

    %% Zu speichernde Variable zusammenstellen
    for i = 1:size(klick_list,2); 
      % Einlesen des Variablennamens und zuweisen seines Inhaltes
      VarName = var_list(klick_list(1,i)).name;

      try
        VarValue = evalin('base',VarName);
      catch
        h=warndlg(['Variable ' VarName ' could not be found and will be ignored henceforth.'],'Scan workspace'); 
        set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
        set(h,'Tag','helpdlgbox','HandleVisibility','On');
        continue;
      end

      % Zeilenvektor in Spaltenvektor transformieren 
      if (size(VarValue,1)) == 1
        VarValue = VarValue';
      end

      % Feststellen der gewünschten Reihen
      SelectRows = klick_list(klick_list(:,i) ~= 0,i);
      SelectRows(1,:) = []; % Lösche Kopf (Enthält Verweis auf den Variablennamen)

      % Wenn Spalten gewählt wurden, dann den Varibleninhalt darauf begrenzen
      if ~isempty(SelectRows)
         VarValue = VarValue(:,SelectRows);
      end

      % Unterschiedliche Länge -> Auffüllen mit "NaN"
      if (size(Matrix,1) ~= size(VarValue,1)) & (size(Matrix,1) ~= 0)
        % Matrix zu lang
        if (size(Matrix,1) > size(VarValue,1))  
           VarValue((size(VarValue,1)+1):size(Matrix,1),:) = NaN;
        % Spaltenvektor zu lang
        else
           Matrix((size(Matrix,1)+1):size(VarValue,1),:) = NaN;
        end
      end

      Matrix = [Matrix VarValue]; 
    end

    %% Speichern
    eval('SaveMatrix = Matrix;');
    save([Pathname Filename], 'SaveMatrix', '-ascii');
    setptr(gcf,'arrow')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case('add_programme')

    setptr(gcf,'watch')
    prompt={'Name of the programme to be added:','Minimal number of parameters (e.g. 1):','Maximal number of parameters (max. 4):'};
    def={'','',''};
    set(0,'DefaultUIControlBackgroundColor',props.msgbox.BackgroundColor);
    answer = inputdlg(prompt,'Add a Matlab programme',1,def);
    set(0,props.root);
    if ~isempty(answer) 
       if isempty(which(answer{1})) & ~isempty(answer{1})
          h=errordlg(['Programme not found.',char(10),'Please ensure, that the programme is inside the Matlab search path.'],'Add a Matlab programme');
          set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
          waitfor(h)
          setptr(gcf,'arrow')
          return
       elseif isempty(answer{1})
          setptr(gcf,'arrow')
          return
       else
          if answer{2}=='', answer{2}='0'; end
          if answer{3}=='', answer{3}='0'; end
          if str2double(answer{3})>4, answer{3}='4'; 
             h=warndlg(['In this version, the number of parameters is',char(10),...
                        'limited to 4. Your input was set to this value.'],'Add a Matlab programme'); 
             set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
             set(h,'Tag','helpdlgbox','HandleVisibility','On');
             waitfor(h)
          end
          if str2double(answer{2})>4, answer{2}='4'; 
             h=warndlg(['In this version, the number of parameters is',char(10),...
                        'limited to 4. Your input was set to this value.'],'Add a Matlab programme'); 
             set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
             set(h,'Tag','helpdlgbox','HandleVisibility','On');
             waitfor(h)
          end
          if str2double(answer{2})>str2double(answer{3}), answer{2}=answer{3}; 
             h=warndlg(['The maximal number of parameters can not be smaller',char(10),...
                        'than the minimal number of parameters. The minimal',char(10),...
                        'number of parameters was set to ',char(answer{3}),'.'],'Add a Matlab programme'); 
             set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
             set(h,'Tag','helpdlgbox','HandleVisibility','On');
             waitfor(h)
          end
          flag=0;
          for i=1:length(funktion),
             if strcmpi(funktion(i).name,answer{1})
                h=errordlg('Programme already installed.','Add a Matlab programme');
                set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
                waitfor(h)
                flag=1;
                break
             end
          end
          if flag==0
            mguircfile=which('mgui.rc');
            fid=fopen(mguircfile,'a');
            if fid>1;
              fprintf(fid,'%s\n',[answer{1},' ',answer{2},' ',answer{3}]);
            else
                h=errordlg('Could not open the mgui.rc file.','Add a Matlab programme');
                set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
                waitfor(h)
            end
            fclose(fid);
            funktion(end+1).name = answer(1);
            funktion(end).min = str2double(answer{2});
            funktion(end).max = str2double(answer{3});
            funktion(end).akt = str2double(answer{2}); 
            funktion(end).klicklist{4}=[];
            funktion(end).klickstring{4}=[];
            set(hfs,'UserData',funktion,'String', [funktion.name]);
%            set(hfs,'Value',length(funktion));
%            mgui mod2;
            set(0,'DefaultUIControlBackgroundColor',props.msgbox.BackgroundColor);
            answer2=questdlg([upper(char(answer{1})),' was successfully installed and',char(10),... 
                'is now added to the function listbox.'],'Add a Matlab programme',...
                'OK','Show mgui.rc','OK'); 
            set(0,props.root);
            if strcmp(answer2,'Show mgui.rc')
               if isunix
                 editor='vi';
               else
                 editor='notepad';
               end
               builtinEd=0;
               
               if isempty(javachk('mwt', 'The MATLAB Editor'));
                 if com.mathworks.ide.editor.EditorOptions.getBuiltinEditor == 0,
                    editor = char(com.mathworks.ide.editor.EditorOptions.getOtherEditor);
                 else
                    com.mathworks.ide.editor.EditorApplication.openDocument(mguircfile);
                    builtinEd=1;
                 end;
               end
               
               if isunix & builtinEd==0
                  if strcmp(editor,'vi') == 1
                     editor = 'xterm -e vi';
                  end
                  eval(['!' editor ' "' mguircfile '" &'])
               elseif ~isunix & builtinEd==0
                 eval(['!"' editor '" "' mguircfile '" &'])
               end
                              
            end
          end
       end
    end 
    setptr(gcf,'arrow')
end

set(0,props.root);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lokal_ws=check_workspace
    %% Aktualisieren der Variablen
    % Einlesen aller vorhandenen VariablenNamen und ihre Feldgrösse
    lokal_ws_all = evalin('base','whos');
    hMF = findobj('Tag','MainFrame');
    hfs = findobj('Tag','function_select');
    if ~(isempty(hMF)), props = get(hMF,'UserData'); else error('Properties not found.'), end
    select = get(hfs,'Value');
    funktion = get(hfs,'UserData');
    FN_Name = funktion(select).name;

    if length(lokal_ws_all) > 0
      j = 1;
      for i=1:length(lokal_ws_all)
        if strcmp(lokal_ws_all(i).class,'double') | strcmp(lokal_ws_all(i).class,'single') | strcmp(lokal_ws_all(i).class,'int8') | strcmp(lokal_ws_all(i).class,'uint8') | strcmp(lokal_ws_all(i).class,'int16') | strcmp(lokal_ws_all(i).class,'uint16') | strcmp(lokal_ws_all(i).class,'int32') | strcmp(lokal_ws_all(i).class,'uint32')
          Ausgabe(j,1) = cellstr([lokal_ws_all(i).name, ' [', num2str(lokal_ws_all(i).size(1)),'x',num2str(lokal_ws_all(i).size(2)), ']' ]);
          lokal_ws(j) = lokal_ws_all(i);
          j = j+1;
        end
      end 
    % Der Workspace ist leer 
    end
    if ~exist('lokal_ws','var') | isempty(lokal_ws)
      if ~strcmp(FN_Name,'load')
         h=errordlg('No numeric variables in the workspace.','Scan workspace');
         set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
         set(h,'Tag','helpdlgbox','HandleVisibility','On');
      end
      Ausgabe(1,1) = cellstr('Empty');
      lokal_ws(1).name = [];
    end
    neu_var={lokal_ws.name};

 for i=1:4;
   hparm = findobj('Tag',['button_parm' num2str(i)]);
   hlbox = findobj('Tag',['listbox_parm' num2str(i)]);
   for l=5:length(funktion),
      temp=funktion(l).klicklist;
      if ~isempty(temp), klick_list=temp{i}; else klick_list = []; end
      temp=get(hlbox,'UserData');
      if ~isempty(temp) & ~isempty(klick_list)
        old_var={temp.name};
        old_klicklist=klick_list(1,:);
        for k=1:length(old_klicklist)
          if length(old_klicklist)<k, klick_list(1,k)=NaN; break, end
          if length(old_var)<old_klicklist(k), klick_list(1,k)=NaN; break, end
          temp=find(strcmp(neu_var,old_var(old_klicklist(k))));
          if ~isempty(temp)
            klick_list(1,k)=temp;
          else
            klick_list(1,k)=NaN;
          end
        end
        klick_list(:,isnan(klick_list(1,:)))=[];
        klick_string=make_klick_string(klick_list,lokal_ws);
        funktion(l).klicklist{i}=klick_list;
        funktion(l).klickstring{i}=klick_string;
        if select==l,
           set(hlbox,'Value',klick_list(1,:));
           set(hparm,'String',klick_string);
           set(hparm,'ToolTip',klick_string);
         end
      end
   end
   set(hfs,'UserData',funktion)
   set(hlbox,'String',Ausgabe,'UserData',lokal_ws);
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function klick_string=make_klick_string(klick_list,var_list)

    klick_string='';

    if ~isempty(klick_list)

      for i = 1:size(klick_list,2)
        klick_string = strcat(klick_string,var_list(klick_list(1,i)).name);
        selected = klick_list(klick_list(:,i) ~= 0,i);

        if size(selected,1) > 1 % Wenn Einzelwerte gefunden wurden
          klick_string = strcat(klick_string, '('); 

          for j = 2:size(selected,1)
            klick_string = strcat(klick_string,num2str(selected(j)));
            klick_string = strcat(klick_string, ','); 
          end
          klick_string(end) = []; 

          klick_string = strcat(klick_string, ')'); 
        end 
 
        klick_string = strcat(klick_string, ','); 
      end
      klick_string(end) = []; 

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check_apply

hfs    = findobj('Tag','function_select');
select = get(hfs,'Value');
funktion = get(hfs,'UserData');

    %% Alle "Klick_listen" durchsuchen, ob "Apply" ausgeführt werden darf
    FN_Akt   = funktion(select).akt;

    trigger = 0;
    % Alle "Klick_list"en durchsuchen
    for i = 1 : FN_Akt
      hbcp = findobj('Tag',['button_close_parm' num2str(i)]);
      temp=funktion(select).klicklist;
      if ~isempty(temp), klick_list=temp{i}; else klick_list=[]; end
      if isempty(klick_list), klick_list=get(hbcp,'UserData'); end
      if isempty(klick_list)
        trigger = 1;
      end
    end
    % Wenn mindestens eine leere dabei, dann die Ausführung verhindern
    if trigger
      set(findobj('Tag','apply'),'Enable','off');
    else
      set(findobj('Tag','apply'),'Enable','on');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ENDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
