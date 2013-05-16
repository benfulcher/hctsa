function xout=jrqa(varargin)
%JRQA   Computes and plots the JRQA measures.
%    Y=JRQA(X [,Y] [,param1,param2,...]) 
%    Recurrence quantification analysis of the joint recurrence
%    plot of the vectors X and Y. 
%
%    The input vectors can be multi-column vectors, where
%    each column will be used as a component of the 
%    phase-space vector. However, if the first column is
%    monotonically increasing, it will be used as an
%    time scale for plotting.
%
%    Y=JRQA(X,M,T,E,W,WS,LMIN,VMIN,TW) computes the 
%    recurrence quantification analysis of the recurrence
%    plot of X by using the dimension M, delay T, the
%    size of neighbourhood E, the window size W and 
%    a window shifting value of WS. LMIN and VMIN 
%    specify the minimal length of diagonal and vertical 
%    line structures (default is 2) and TW specifies the
%    Theiler window (default is 1).
%
%    JRQA(...) without any output arguments opens a
%    GUI for interactively control the JRQA. If an
%    output is specified with using the option 'gui',
%    then the output will contain the figure handle.
% 
%    Parameters: 
%    Dimension M, delay T, the size of
%    neighbourhood E, the window size W and the shift
%    value WS are the first five numbers after the data 
%    series; if W is empty, the whole plot will be calculated.
%    The next two optional numeric parameters LMIN and VMIN
%    specify the minimal length of line structures.
%
%    As the last numeric parameter, the size of the Theiler 
%    window TW can be specified (default is 1). This window 
%    excludes the recurrence points parallel to the main 
%    diagonal from the analysis. The application of the
%    Theiler window is useful only for recurrence plots. In
%    joint recurrence plots, the size of the Theiler window will
%    be set automatically to zero.
%
%    Further parameters can be used to switch between various 
%    methods of finding the neighbours of the phasespace 
%    trajectory, to suppress the normalization of the data 
%    and to suppress the GUI (useful in order to use this 
%    programme by other programmes). The minimal length of 
%    diagonal and vertical structures can be setted only in 
%    the GUI.
%
%    Methods of finding the neighbours.
%      maxnorm     - Maximum norm.
%      euclidean   - Euclidean norm.
%      minnorm     - Minimum norm.
%      nrmnorm     - Euclidean norm between normalized vectors
%                    (all vectors have the length one).
%      rr          - Maximum norm, fixed recurrence rate.
%      fan         - Fixed amount of nearest neighbours.
%      inter       - Interdependent neighbours.
%      omatrix     - Order matrix.
%      opattern    - Order patterns recurrence plot.
%
%    Normalization of the data series.
%      normalize   - Normalization of the data.
%      nonormalize - No normalization of the data.
%
%    Suppressing the GUI.
%      gui         - Creates the GUI.
%      nogui       - Suppresses the GUI.
%      silent      - Suppresses all output.
%
%    Parameters not needed to be specified.
%
%    The window of length w is applied on the data and not on the RP, 
%    i.e. the RP will have smaller size than the window, thus w-(m-1)*tau. 
%    If we consider the data window at time i ... i+w, the corresponding RQA 
%    measures are assigned to time i. Therefore, if you see a beginning 
%    of a transition in the plot of the RQA measures at time i, this 
%    transition will probably happen at time i+w-(m-1)*tau. 
%
%    Output:
%      Y(:, 1) = RR     (recurrence rate)
%      Y(:, 2) = DET    (determinism)
%      Y(:, 3) = <L>    (mean diagonal line length)
%      Y(:, 4) = Lmax   (maximal diagonal line length)
%      Y(:, 5) = ENTR   (entropy of the diagonal line lengths)
%      Y(:, 6) = LAM    (laminarity)
%      Y(:, 7) = TT     (trapping time)
%      Y(:, 8) = Vmax   (maximal vertical line length)
%      Y(:, 9) = T1     (recurrence time of 1st type)
%      Y(:,10) = T2     (recurrence time of 2nd type)
%
%    Examples: N = 500; w = 40; ws = 10;
%              b = .4; a = .6; mu = .8:-0.7/N:.1;
%              % two mutually coupled logistic maps
%              for i = 2:N,
%                  a(i) = 3.6 * a(i-1) * (1 - a(i-1)); 
%                  b(i) = 4 * b(i-1) * (1 - b(i-1)) - mu(i)*a(i); 
%              end
%              plot(b)
%              % coupling is obtained by higher RR and DET values
%              jrqa(a,b,1,1,.2,w,ws);
%
%      
%    See also CRQA, JRP and CRP.
%
%    References: 
%    Trulla, L. L., Giuliani, A., Zbilut, J. P., Webber Jr., C. L.: 
%    Recurrence quantification analysis of the logistic equation with 
%    transients, Phys. Lett. A, 223, 1996.
%
%    Marwan, N., Wessel, N., Meyerfeldt, U., Schirdewan, A., Kurths, J.: 
%    Recurrence Plot Based Measures of Complexity and its Application to 
%    Heart Rate Variability Data, Phys. Rev. E, 66(2), 2002.
%
%    Romano, M., Thiel, M., Kurths, J., von Bloh, W.: 
%    Multivariate Recurrence Plots, Phys. Lett. A , 330, 2004.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 1998-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:32:57 $
% $Revision: 2.7 $
%
% $Log: jrqa.m,v $
% Revision 2.7  2009/03/24 08:32:57  marwan
% copyright address changed
%
% Revision 2.6  2007/07/18 17:18:44  marwan
% integer values in the arguments supported
%
% Revision 2.5  2007/05/15 17:33:13  marwan
% new neighbourhood criterion: fixed RR
%
% Revision 2.4  2006/02/14 11:45:49  marwan
% *** empty log message ***
%
% Revision 2.3  2005/11/28 10:16:35  marwan
% && and || changed to & and |
% (seems to cause problems in Matlab 12.1)
%
% Revision 2.2  2005/04/15 09:02:32  marwan
% minor bugfix in plugin section
%
% Revision 2.1  2005/04/08 09:54:19  marwan
% plugin added
%
% Revision 1.1  2005/03/16 16:22:29  marwan
% support for joint recurrence plots added
%
%
%
% This program is part of the new generation XXII series.
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% programme properties

global errcode props


init_properties

lmin=2;
vmin=2;
nonorm=1;
theiler_window=1;
hw=-1;
xscale = [];
yscale = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check the input

error(nargchk(1,13,nargin));
if nargout>1, error('Too many output arguments'), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splash the GPL

splash_gpl('crp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check and read the input

%try 
errcode=1;

   set(0,'ShowHidden','on')
   varargin{14}=[];
   i_double=find(cellfun('isclass',varargin,'double'));
   i_char=find(cellfun('isclass',varargin,'char'));
   t=1;
   w=[];wstep=0; method='max'; method_n=1;
   time_scale_flag=1; time_scale_flag_x=1;
   nogui=0;
      check_meth={'ma','eu','mi','nr','rr','fa','in','om','op','di'}; 	% maxnorm, euclidean, nrmnorm,  fan, distance
      check_gui={'gui','nog','sil'};                % gui, nogui, silent
      check_norm={'non','nor'};                        % nonormalize, normalize
   
   if nargin & ischar(varargin{1})
     set(0,'ShowHidden','on');h=findobj('Tag','Msgbox_Check Input');delete(h)
     action=varargin{1};
     h=findobj('Tag','crqa_Fig');
     temp=get(h(1),'UserData');
     x=temp{1}; y=temp{2};
     h=findobj('Tag','crqa_axes_Data','Parent',h(1)); % get handle of plot-axes object
     h_crqa_axes_Data = h(1);
     h=findobj('Type','line','Parent',h(1)); % get handle of line-object
     xscale=get(h(1),'xdata');
     h=findobj('Tag','crqa_m');
     m=str2num(get(h(1),'String'));
     h=findobj('Tag','crqa_maxLag');
     t=str2num(get(h(1),'String'));
     h=findobj('Tag','crqa_eps');
     e=str2num(get(h(1),'String'));
     h=findobj('Tag','crqa_method');
     method={'Maximum Norm','Euclidean Norm','Minimum Norm','Normalized Norm','Maximum Norm, fixed RR','FAN','Interdependent','Order Matrix','Order Pattern','Distance Plot'};
     method=method{get(h(1),'Value')};
     nonorm=get(h(1),'UserData');
     h=findobj('Tag','crqa_lmin');
     lmin=str2num(get(h(1),'String'));
     h=findobj('Tag','crqa_vmin');
     vmin=str2num(get(h(1),'String'));
     h=findobj('Tag','crqa_theiler');
     theiler_window=str2num(get(h(1),'String'));
     h=findobj('Tag','crqa_w');
     w=str2num(get(h(1),'String'));
     h=findobj('Tag','crqa_ws');
     wstep=str2num(get(h(1),'String'));

   elseif nargin & isnumeric(varargin{1})

     % check the text input parameters for method, gui 
      temp_meth=0;
      temp_norm=0;
      temp_gui=0;
      if ~isempty(i_char)
         for i=1:length(i_char), 
            varargin{i_char(i)}(4)='0';
            temp_gui=temp_gui+strcmpi(varargin{i_char(i)}(1:3),check_gui'); 
            temp_norm=temp_norm+strcmpi(varargin{i_char(i)}(1:3),check_norm'); 
            temp_meth=temp_meth+strcmpi(varargin{i_char(i)}(1:2),check_meth'); 
         end
         method_n=min(find(temp_meth));
         nogui=min(find(temp_gui))-1;
         nonorm=min(find(temp_norm))-1;
         for i=1:length(i_char); temp2(i,:)=varargin{i_char(i)}(1:3); end
         i_char(strmatch(check_gui(find(temp_gui)),temp2))=[];
         if isempty(nonorm), nonorm=1; end
         if nonorm>1, nonorm=1; end
         if isempty(nogui), nogui=0; end
         if isempty(method_n), method_n=1; end
         if nogui>2, nogui=1; end
         if method_n>length(check_meth), method0=length(check_meth); end
         method=check_meth{method_n};
      else
         nogui=0;
         if nargout
            nogui=1;
            action='compute'; 
         end
      end
      if nogui==0
        action='init';
      else
        action='compute'; 
      end
      
      % get the parameters for creating RP
      if max(size(varargin{1}))<=3
         disp('Error using ==> jrqa')
         disp('To less values in data X.')
         return
      end
      x=double(varargin{1});
      if isempty(varargin{2}) | ~isnumeric(varargin{2}), y=x; else
      y=double(varargin{2}); end
%      if sum(double(diff(x(:,1))<=0)), time_scale_flag=0; end
   
      if (isnumeric(varargin{2}) & max(size(varargin{2}))==1) | ~isnumeric(varargin{2})
        y=x;
        if ~isempty(varargin{i_double(2)}), m=varargin{i_double(2)}(1); else m=1; end
        if ~isempty(varargin{i_double(3)}), t=varargin{i_double(3)}(1); else t=1; end
        if ~isempty(varargin{i_double(4)}), e=varargin{i_double(4)}(1); else e=.1; end
        if ~isempty(varargin{i_double(5)}), w=varargin{i_double(5)}(1); else w=varargin{i_double(5)}; end
        if ~isempty(varargin{i_double(6)}), wstep=varargin{i_double(6)}(1); else wstep=1; end
        if ~isempty(varargin{i_double(7)}), lmin=varargin{i_double(7)}(1); end
        if ~isempty(varargin{i_double(8)}), vmin=varargin{i_double(8)}(1); end
        if ~isempty(varargin{i_double(9)}), theiler_window=varargin{i_double(9)}(1); end
      else
        if ~isempty(varargin{i_double(3)}), m=varargin{i_double(3)}(1); else m=1; end
        if ~isempty(varargin{i_double(4)}), t=varargin{i_double(4)}(1); else t=1; end
        if ~isempty(varargin{i_double(5)}), e=varargin{i_double(5)}(1); else e=.1; end
        if ~isempty(varargin{i_double(6)}), w=varargin{i_double(6)}(1); else w=varargin{i_double(6)}; end
        if ~isempty(varargin{i_double(7)}), wstep=varargin{i_double(7)}(1); else wstep=1; end
        if ~isempty(varargin{i_double(8)}), lmin=varargin{i_double(8)}(1); end
        if ~isempty(varargin{i_double(9)}), vmin=varargin{i_double(9)}(1); end
        if ~isempty(varargin{i_double(10)}), theiler_window=varargin{i_double(10)}(1); end
      end
    else
      disp('Error using ==> jrqa')
      disp('No valid arguments.')
      return
    end

   
   if max(size(x))~=max(size(y)),
        if ~nogui, errordlg('Data must have the same length.','Check Data'), waitforbuttonpress, return, else error('Data must have the same length.'), end
   end
   
    Nx=length(x); Ny=length(y);
    if size(x,1)<size(x,2), x=x'; end
    if size(y,1)<size(y,2), y=y'; end
   
   if size(x,2)>=2
      xscale=x(:,1); 
      if ~isempty(find(diff(xscale)<0)) % multi-column data vector, each column used as vector component
          time_scale_flag=0; 
          time_scale_flag_x=0; 
          xscale=(1:length(x))'; 
      end
   else
      xscale=(1:length(x))'; 
      time_scale_flag=0; 
      time_scale_flag_x=0; 
   end

   if size(y,2)>=2
      yscale=y(:,1); 
      time_scale_flag=1; 
      if ~isempty(find(diff(yscale)<0))
          time_scale_flag=0; 
          yscale=(1:length(y))';
      end
      if time_scale_flag & ~time_scale_flag_x % if time-scale given in y, but not in x -> error
          if ~nogui
             errordlg(['A time-scale for the second data is series given, but not for the first!',10,'(The time-scale has to be inlcuded as the first colummn of the first data vector.)'],'Check Data')
             waitforbuttonpress
             return
          else
             error(['A time-scale for the second data series is given, but not for the first!',10,'(The time-scale has to be inlcuded as the first colummn of the first data vector.)'])
          end
      end
   else
       yscale=(1:length(y))';
       time_scale_flag=0; 
   end
   if time_scale_flag_x & ~time_scale_flag % if time-scale given in x, but not in y
       if length(x) ~= length(y)
          if ~nogui
             errordlg(['If you are using the time-scale given by the first vector also for',10,'the second vector, both vectors should have the same size!'],'Check Data')
             waitforbuttonpress
             return
          else
             error(['If you are using the time-scale given by the first vector also for',10,'the second vector, both vectors should have the same size!'])
          end
       end
       y = [xscale, y];
       yscale = xscale;
       time_scale_flag=1; 
   end
   if e<0, 
        e=1; 
        if ~nogui
           warndlg('The threshold size E can not be negative and is now set to 1.','Check Data')
           waitforbuttonpress
           h=findobj('Tag','crqa_eps');
           if ~isempty(h), set(h(1),'String',num2str(e)), end
        else 
           disp('The threshold size E can not be negative and is now set to 1.'), 
        end
      end
      if t<1, 
        t=1; 
        if ~nogui
           warndlg('The delay T can not be smaller than one and is now set to 1.','Check Data')
           waitforbuttonpress
           h=findobj('Tag','crqa_maxLag');
           if ~isempty(h), set(h(1),'String',num2str(t)), end
        else
           disp('The delay T can not be smaller than one and is now set to 1.')
        end
      end
      if isempty(w), w=Nx; wstep=1; end
      if w < 5+(m-1)*t, 
        w=5+(m-1)*t;
        if ~nogui, warndlg('The window size W exceeds the valid range.','Check Data')
           waitforbuttonpress
           h=findobj('Tag','crqa_w');
           if ~isempty(h), set(h(1),'String',num2str(w)), end
        else, disp('The window size W exceeds the valid range.'), end
      end
      if w>Nx, 
        w=Nx; wstep=1;; 
        if ~nogui, warndlg('The window size W exceeds the valid range.','Check Data')
           waitforbuttonpress
           h=findobj('Tag','crqa_w');
           if ~isempty(h), set(h(1),'String',num2str(w)), end
        else, disp('The window size W exceeds the valid range.'), end
      end
      if wstep<1 | wstep>Nx/3, 
        wstep=2; 
        if ~nogui, warndlg('The window shifting value WS exceeds the valid range.','Check Data')
           waitforbuttonpress
           h=findobj('Tag','crqa_ws');
           if ~isempty(h), set(h(1),'String',num2str(wstep)), end
        else, disp('The window shifting value WS exceeds the valid range.'), end
      end
      if vmin<1 | vmin>Nx, 
        vmin=2; 
        if ~nogui, warndlg('The minimal length for vertical lines is not valid.','Check Data')
           waitforbuttonpress
           h=findobj('Tag','crqa_vmin');
           if ~isempty(h), set(h(1),'String',num2str(vmin)), end
        else, disp('The minimal length for vertical lines is not valid.'), end
      end
      if lmin<1 | lmin>Nx, 
        lmin=2; 
        if ~nogui, warndlg('The minimal length for diagonal lines is not valid.','Check Input')
           waitforbuttonpress
           h=findobj('Tag','crqa_lmin');
           if ~isempty(h), set(h(1),'String',num2str(lmin)), end
        else, disp('The minimal length for diagonal lines is not valid.'), end
      end
      if theiler_window<0 | theiler_window>Nx, 
        theiler_window=1; 
        if ~nogui, warndlg('The value for the Theiler window is not valid.','Check Input')
           waitforbuttonpress
           h=findobj('Tag','crqa_theiler');
           if ~isempty(h), set(h(1),'String',num2str(theiler_window)), end
        else, disp('The value for the Theiler window is not valid.'), end
      end
      if x~=y, theiler_window=0; end
      t=round(t); m=round(m); w=round(w); wstep=round(wstep); vmin=round(vmin); lmin=round(lmin); theiler_window=round(theiler_window);

switch(action)

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% create gui

case 'init'

  errcode=3;
       
  oldunit=get(0,'unit');
  set(0,'Unit','char')
  scr=get(0,'ScreenSize'); 
  set(0,'Unit',oldunit)
  
  h8=figure(props.window,...                                % Plot Figure
            'Tag','crqa_Fig',...
            'MenuBar','Figure',...
            'Position',[(scr(3)-150)/2 scr(4)-50 150.0000 44],...
            'PaperType','a4',...
            'PaperPosition',[0.25 0.25 7.7677 11.193],...
            'PaperOrientation','portrait',...
            'UserData',{x,y,time_scale_flag},...
              'Name','Joint Recurrence Quantification Analysis');
  
  set(0,'showhidden','on')
  h=findobj('Label','&Help','Type','uimenu');
  if isempty(h)
    h=uimenu('Label','&Help');
    h2=uimenu('Parent',h(1),'Label','&Help Joint Recurrence Quantification Analysis','Callback','helpwin jrqa');
  else
    h1=flipud(get(h(1),'Children'));
    set(h1(1),'Separator','on')
    h2=uimenu('Parent',h(1),'Label','&Help Joint Recurrence Quantification Analysis','Callback','helpwin jrqa');
    copyobj(h1,h(1))
    delete(h1)
  end
  set(0,'showhidden','off')


  h=axes(props.axes,...
            'Position',[89+30 24.8+14.5 6.8 3.5]);    
  logo=load('logo');
  h2=imagesc([logo.logo fliplr(logo.logo)]);
  set(h2,'Tag','uniLogo')
  set(h,props.logo,'Tag','axes_logo')
  h=uicontrol(props.text,...
            'Tag','text_logo',...
            'String','Uni Potsdam',...
            'Position',[97+30 24.2143+14.5 22  3.5714]);    
  h2=textwrap(h,{[char(169),' AGNLD'],'University of Potsdam','1998-2007'});
  set(h,'String',h2)


%%%%%%%%%%% plots

  axes_height = 5.8;
  axes_base = 3;
  axes_hoffset = 2.5;
  h=axes(props.axes,...
            'Tag','crqa_axes_Data',...
                'Box','On',...
            'Position',[11.2017  axes_base+4*(axes_height+axes_hoffset)    28.1785+15    axes_height]);    

  if time_scale_flag
      plot(xscale,x(:,2:end),'color',props.line.Color)
      if ~all(x(:)==y(:)) hold on; plot(yscale,y(:,2:end),'r'), end
  else
      plot(xscale,x(:,1),'color',props.line.Color)
      if ~all(x(:)==y(:)) hold on; plot(yscale,y(:,1),'r'), end
  end
  set(h,'Tag','crqa_axes_Data','color',props.axes.Color)
  ylabel('Data')

  h=axes(props.axes,...
            'Tag','crqa_axes_Var',...
            'Box','On',...
            'Position',[49.8023+15   axes_base+4*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  if time_scale_flag
     x_var = winplot(x(:,2:end),w,wstep,2);
  else
     x_var = winplot(x(:,1),w,wstep,2);
  end
  if size(x_var,1) > 2
    h2=stairs(xscale(round(x_var(:,1))),x_var(:,2));
    set(h2,'color',props.line.Color)
  else
    cla
    text(0.5,0.5,sprintf('%6.4f',x_var(1,2)),'FontWeight','bold','HorizontalAlign','Center')
  end
  if ~all(x(:)==y(:)) 
    hold on; 
    if time_scale_flag
        y_var = winplot(y(:,2:end),w,wstep,2); 
    else
        y_var = winplot(y(:,1),w,wstep,2); 
    end
    if size(y_var,1) > 2
      stairs(yscale(round(y_var(:,1))),y_var(:,2),'r')
    else
      cla
      text(0.5,0.5,sprintf('%6.5f, %6.5f',x_var(2,2),y_var(1,2)),'FontWeight','bold','HorizontalAlign','Center')
    end
    h3=axes('Units','Char',...
              'Pos',get(h,'Pos'),...
              'XLim',get(h,'XLim'),...
              'YLim',get(h,'YLim'),...
              'YAxisLocation','right',...
              'Color','none',...
              'Tag','crqa_axes_CoVar',...
              'visible','off');
  end
  set(gcf,'CurrentAxes',h);
  ylabel('Variance')
  set(h,'Tag','crqa_axes_Var','color',props.axes.Color)


  h=axes(props.axes,...
            'Tag','crqa_axes_RR',...
            'Box','On',...
            'Position',[11.2017   axes_base+3*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('RR')

  h=axes(props.axes,...
            'Tag','crqa_axes_DET',...
            'Box','On',...
            'Position',[49.8023+15   axes_base+3*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('DET')

  h=axes(props.axes,...
            'Tag','crqa_axes_L',...
            'Box','On',...
            'Position',[11.2017   axes_base+2*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('L')

  h=axes(props.axes,...
            'Tag','crqa_axes_ENTR',...
            'Box','On',...
            'Position',[49.8023+15   axes_base+2*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('ENTR')

  h=axes(props.axes,...
            'Tag','crqa_axes_LAM',...
            'Box','On',...
            'Position',[11.2017    axes_base+1*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('LAM')

  h=axes(props.axes,...
            'Tag','crqa_axes_TT',...
            'Box','On',...
            'Position',[49.8023+15    axes_base+1*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('TT')

  h=axes(props.axes,...
            'Tag','crqa_axes_T1',...
            'Box','On',...
            'Position',[11.2017    axes_base+0*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('T_1')

  h=axes(props.axes,...
            'Tag','crqa_axes_T2',...
            'Box','On',...
            'Position',[49.8023+15    axes_base+0*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('T_2')

%%%%%%%%%%% embedding
  h=uicontrol(props.frame,...
            'Tag','frame',...
            'Position',[86+30 29.+2.8 29 5.7]);    
 
  h=uicontrol(props.text,...
            'Tag','text',...
            'Fontangle','italic',...
            'String','Embedding parameters',...
            'Position',[87+30 34.2+1.6 16.8333  1.5000]);    
 
  h=uicontrol(props.text,...
            'Tag','text',...
            'String','Dimension:',...
            'Position',[89+30 32+2.2 16.8333  1.5000]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_m',...
            'String',num2str(m),...
                'ToolTip','Select the embedding dimension.',...
            'Position',[104+30 32+.2+2.2 7  1.5000]);    

  h=uicontrol(props.text,...
            'Tag','text',...
            'String','Delay:',...
            'Position',[89+30 30+2.32 16.8333  1.5000]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_maxLag',...
            'String',num2str(t),...
            'ToolTip','Insert the embedding delay time.',...
            'Position',[104+30 30+.2+2.32 7  1.5000]);    

%%%%%%%%%%% neigbourhood
  h=uicontrol(props.frame,...
            'Tag','frame',...
            'Position',[86+30 21.9+2.5 29 6.7]);    
 
  h=uicontrol(props.text,...
            'Tag','text',...
            'Fontangle','italic',...
            'String','Neighbourhood',...
            'Position',[87+30 26.8+2.5 23  1.5000]);    
 
  h=uicontrol(props.popup,...
            'Tag','text',...
            'Tag','crqa_method',...
            'UserData',nonorm,...
            'Value',method_n,...
            'String','Maximum Norm|Euclidean Norm|Minimum Norm|Normalized Norm|Fixed RR|Fixed Amount|Interdependent|Order Matrix|Order Patterns|Distance Plot',...
            'Position',[89+30 24.6+.2+2.5 22  1.7]);    

  h=uicontrol(props.text,...
            'Tag','text',...
            'String','Threshold:',...
            'Position',[89+30 22.6+2.5 16.8333  1.5]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_eps',...
            'String',num2str(e),...
            'ToolTip','Insert the size of neighbourhood.',...
            'Position',[104+30 22.6+.2+2.5 7  1.5000]);    

%%%%%%%%%%% jrqa parameters
  h=uicontrol(props.frame,...
            'Tag','frame',...
            'Position',[86+30 11+.5 29 12.2]);    
 
  h=uicontrol(props.text,...
            'Tag','text',...
            'Fontangle','italic',...
            'String','JRQA parameters',...
            'Position',[87+30 19.4+2.5 23  1.5000]);    

  h=uicontrol(props.text,...
            'Tag','text',...
            'String','min. Diagonal:',...
            'Position',[89+30 17.6+2.5 16.8333  1.5]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_lmin',...
            'String',num2str(lmin),...
            'ToolTip','Insert the minimal length of a diagonal line.',...
            'Position',[104+30 17.6+.2+2.5 7  1.5000]);    

  h=uicontrol(props.text,...
            'Tag','text',...
            'String','min. Vertical:',...
            'Position',[89+30 15.6+2.5 16.8333  1.5]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_vmin',...
            'String',num2str(vmin),...
            'ToolTip','Insert the minimal length of a vertical line.',...
            'Position',[104+30 15.6+.2+2.5 7  1.5000]);    

  h=uicontrol(props.text,...
            'Tag','text',...
            'String','Theiler wind.:',...
            'Position',[89+30 15.6+.5 16.8333  1.5]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_theiler',...
            'String',num2str(theiler_window),...
            'ToolTip','Insert the size for the Theiler window.',...
            'Position',[104+30 15.6+.2+.5 7  1.5000]);    
  if x~=y; set(h,'enable','off'); end

  h=uicontrol(props.text,...
            'Tag','text',...
            'String','Window size:',...
            'Position',[89+30 13.6+.5 16.8333  1.5]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_w',...
            'String',num2str(w),...
            'ToolTip','Insert the size of the sliding window.',...
            'Position',[104+30 13.6+.2+.5 7  1.5000]);    

  h=uicontrol(props.text,...
            'Tag','text',...
            'String','Window step:',...
            'Position',[89+30 11.6+.5 16.8333  1.5]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_ws',...
            'String',num2str(wstep),...
            'ToolTip','Insert the step width for sliding the window.',...
            'Position',[104+30 11.6+.2+.5 7  1.5000]);    


%%%%%%%%%%% buttons

  h=uicontrol(props.frame,...
            'Tag','frame',...
            'Position',[86+30 1.3+.5 29 9]);    


 h=uicontrol(props.button,...
                'String','Store',...
                'Tag','crqa_button_store',...
              'Enable','Off',...
                'ToolTip','Stores the JRQA analysis into a variable in the workspace.',...
                'Callback','jrqa store',...
                'Position',[100.5+30  7.4+.5 10.5  2.2143]);

  h=uicontrol(props.button,...
            'Tag','crqa_button_print',...
            'CallBack','jrqa print',...
                'ToolTip','Prints the JRQA window.',...
            'String','Print',...
            'Position',[89+30  7.4+.5 10.5  2.2143]);    

  h=uicontrol(props.button,...
            'Tag','crqa_button_close',...
            'CallBack','jrqa close',...
                'ToolTip','Closes the JRQA window.',...
            'String','Close',...
            'Position',[89+30  2.1+.5 22  2.2143]);    
  
     h=uicontrol(props.button,...
                'String','Apply',...
                'Tag','crqa_button_apply',...
                'ToolTip','Starts the computation.',...
                'Callback','jrqa compute',...
                'Position',[89+30  4.75+.5 22  2.2143]);
  
  set(0,'ShowHidden','on')
  set(h8, 'HandleVis','CallBack')
  tags={'crqa_Fig';'axes_logo';'text_logo';'crqa_theiler';'frame';'text';'crqa_axes_Data';'crqa_axes_Var';'crqa_axes_CoVar';'crqa_axes_RR';'crqa_axes_DET';'crqa_axes_L';'crqa_axes_ENTR';'crqa_axes_LAM';'crqa_axes_TT';'crqa_axes_T1';'crqa_axes_T2';'crqa_m';'crqa_maxLag';'crqa_method';'crqa_eps';'crqa_lmin';'crqa_vmin';'crqa_w';'crqa_ws';'crqa_button_store';'crqa_button_print';'crqa_button_close';'crqa_button_apply'};
  h=[];
  for i=1:length(tags); h=[h; findobj('Tag',tags{i})]; end
  set(h,'Units','Norm')
  if nargout; xout=h8; end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% close windows 

case 'close'
  errcode=101;
  set(0,props.root)
  h=findobj('Tag','crqa_Fig');
  if ~isempty(h), close(h(1)), end 
  clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% store

case 'store'
  errcode=0;
  if ~isempty(findobj('Tag','crqa_button_store'))
    h=findobj('Tag','crqa_button_store');
    h1=findobj('Tag','crqa_button_close');
    if ~isempty(h1), vname_old=get(h1(1),'UserData'); else vname_old=''; end
    if isempty(vname_old), vname_old=''; end
    vname=char(inputdlg('Choose a variable name.','Store output',1,{vname_old}));
    if isempty(vname)
      return
    else
      crqa_values=get(h(1),'UserData');
      assignin('base',vname, [crqa_values])
      warndlg(['JRQA measures have been assigned to the workspace variable ''',vname,'''.'],'Store output');
           waitforbuttonpress
      set(h1(1),'UserData',vname)
    end
  end
  set(0,props.root)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% print
  case 'print'

    errcode=91;
    h=findobj('Tag','crqa_axes_Data','Parent',gcf); h_axes.h(1)=h(1);
    h=findobj('Tag','crqa_axes_Var','Parent',gcf); h_axes.h(2)=h(1);
    h=findobj('Tag','crqa_axes_CoVar','Parent',gcf); h_axes.h(11)=h(1);
    h=findobj('Tag','crqa_axes_RR','Parent',gcf); h_axes.h(3)=h(1);
    h=findobj('Tag','crqa_axes_DET','Parent',gcf); h_axes.h(4)=h(1);
    h=findobj('Tag','crqa_axes_L','Parent',gcf); h_axes.h(5)=h(1);
    h=findobj('Tag','crqa_axes_ENTR','Parent',gcf); h_axes.h(6)=h(1);
    h=findobj('Tag','crqa_axes_LAM','Parent',gcf); h_axes.h(7)=h(1);
    h=findobj('Tag','crqa_axes_TT','Parent',gcf); h_axes.h(8)=h(1);
    h=findobj('Tag','crqa_axes_T1','Parent',gcf); h_axes.h(9)=h(1);
    h=findobj('Tag','crqa_axes_T2','Parent',gcf); h_axes.h(10)=h(1);
    h=findobj('Tag','uniLogo');
    tags={'text_logo';'frame';'text';'crqa_m';'crqa_maxLag';'crqa_method';'crqa_eps';'crqa_lmin';'crqa_vmin';'crqa_theiler';'crqa_w';'crqa_ws';'crqa_button_store';'crqa_button_print';'crqa_button_close';'crqa_button_apply';};
    for i=1:length(tags); h=[h; findobj('Tag',tags{i},'Parent',gcf)]; end
    set(h,'Visible','Off')
    
    set(h_axes.h,'Units','Character');
    h_axes.old_pos=get(h_axes.h,'Position');

    axes_height = .13;
    axes_base = 0.065;
    axes_hoffset = .06;
    
    for i=2:2:10
    set(h_axes.h(i-1),  'Units','normalize','Position',[0.1300    axes_base+(5-i/2)*(axes_height+axes_hoffset)    0.3270    axes_height])
    set(h_axes.h(i),'Units','normalize','Position',[0.5780    axes_base+(5-i/2)*(axes_height+axes_hoffset)    0.3270    axes_height])
    end
    set(h_axes.h(11),  'Units','normalize','Position',[0.5780    axes_base+(5-2/2)*(axes_height+axes_hoffset)    0.3270    axes_height])
    h_dlg=printdlg;
    waitfor(h_dlg)

    for i=1:10, set(h_axes.h(i),'Units','Character','Position',h_axes.old_pos{i}), set(h_axes.h(i),'Units','Norm'),end
    set(h_axes.h(11),'Units','Character','Position',h_axes.old_pos{11}), set(h_axes.h(11),'Units','Norm')
    set(h,'Visible','On')
    set(0,props.root)







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute
case 'compute'

  errcode=11;
  if length(method)>1 & strcmpi(method(1:2),'di')
      disp('Warning: RQA from distance plot not possible!')
      return
  end
  
  if ~nogui
    h_fig=findobj('tag','crqa_Fig');
    setptr(gcf,'watch'), 
    obj=({'text';'crqa_m';'crqa_maxLag';'crqa_method';'crqa_eps';'crqa_lmin';'crqa_vmin';'crqa_theiler';'crqa_w';'crqa_ws';'crqa_button_store';'crqa_button_print';'crqa_button_close'});
    for j=1:length(obj); 
        h=findobj('Tag',obj{j},'Parent',h_fig(1)); 
        if ~isempty(h)
          set(h,'Enable','Off')
        end
    end
    h=findobj('tag','crqa_button_apply');
    set(h(1),'ToolTip','Stops the computation.','String','Stop','Callback','set(0,''ShowHidden'',''on'');h=findobj(''tag'',''crqa_button_apply'');set(h(1),''String'',''Stopped'');set(0,''ShowHidden'',''off'')')
  end

if Nx==w & wstep<2, wstep=1; Nx=w+1; end
if Nx==w, Nx=w+1; end
if nogui~=2, hw=waitbar(0,['0/',num2str(Nx-w)]);set(hw,'Name','Please Wait!');h1=get(hw,'chil');h1=get(h1,'title'); drawnow; end
if strcmpi(method,'Order Pattern') method = 'op'; end
if strcmpi(method,'Order Matrx') method = 'om'; end
if strcmpi(method,'Maximum Norm, fixed RR') method = 'rr'; end


% check if plugin exist and is executable
[plugin_exist, plugin_name, plugin_path] = is_crp_plugin;
if nogui == 1 &plugin_exist & ( method_n < 4 ) & length(x) == length(y)
    disp('(plugin used)')
end

errcode=20;

for i=1:wstep:Nx-w; 
     if ~nogui
       set(0,'ShowHidden','on')
       h=findobj('tag','crqa_button_apply','Parent',h_fig(1));
       if strcmpi(get(h(1),'string'),'stopped')
         Y(i:Nx-w,1:6)=NaN;
         break
       end
     end

     if time_scale_flag
        x_var=var(x(i:i+w-1,2:end)); y_var=var(y(i:i+w-1,2:end)); 
        temp=cov(x(i:i+w-1,2:end),y(i:i+w-1,2:end));
     else
        x_var=var(x(i:i+w-1,1)); y_var=var(y(i:i+w-1,1)); 
        temp=cov(x(i:i+w-1,1),y(i:i+w-1,1));
     end
     xy_var=temp(1,2);

     do_norm = {'non';'nor'};

     % if plugin exist and method is MAX, MIN or EUC
     
     if plugin_exist & ( method_n < 4 ) & length(x) == length(y) 

         errcode=21;
         warning off
         tmp_xdatafile  = tempname;
         tmp_ydatafile  = tempname;
         tmp_rqadatafile = tempname;
         
         
         % normalize data if necessary
         if nonorm
             x = (x - repmat(mean(x),length(x),1)) ./ repmat(std(x),length(x),1);
             y = (y - repmat(mean(y),length(y),1)) ./ repmat(std(y),length(y),1);
         end
         
         x_tmp = x(i:i+w-1,:);
         y_tmp = y(i:i+w-1,:);
         
         % save data in temporary file
         save(tmp_xdatafile,'x_tmp','-ascii','-tabs');
         save(tmp_ydatafile,'y_tmp','-ascii','-tabs');

         % call extern rp programme
         m_str = {'MAX', 'EUC', 'MIN'};

         [status ] = system([plugin_path,filesep,plugin_name,' -m ',num2str(m), ...
                                       ' -t ',num2str(t), ...
                                       ' -e ',num2str(e), ...
                                       ' -n ',m_str{method_n}, ...
                                       ' -w ',num2str(theiler_window), ...
                                       ' -i ',tmp_xdatafile, ...
                                       ' -j ',tmp_ydatafile, ...
                                       ' -o ',tmp_rqadatafile, ...
                                       ' -J ', ...
                                       ' -s']);
         errcode=22;
         % import RQA
         rqa_in = [];
         try
             fid = fopen(tmp_rqadatafile,'r'); % open RQA data file
             while 1
                 dataStr = fgetl(fid); % read RQA data
                 if ~ischar(dataStr), break, end % leave the loop if end of file
                 if isempty(findstr(dataStr, '#')) % neglect comments line (e.g. header)
                     rqa_in = str2num(dataStr); % import RQA measures into local variable rqa
                 end
             end
             fclose(fid); % close RQA data file

             RR = rqa_in(1);
             DET = rqa_in(2);
             LAM = rqa_in(4);
             Lmax = rqa_in(6);
             L = rqa_in(7);
             ENTR = rqa_in(8);
             Vmax = rqa_in(10);
             TT = rqa_in(11);
             t1 = rqa_in(13);
             t2 = rqa_in(14);
             warning off
             delete(tmp_rqadatafile);
             delete(tmp_xdatafile);
             delete(tmp_ydatafile);
             warning on
         catch 
             warning off
             delete(tmp_rqadatafile);
             delete(tmp_xdatafile);
             delete(tmp_ydatafile);
             warning on
         end

         if nogui~=2 & ishandle(h1), set(h1,'str',[num2str(i),'/',num2str(Nx-w)]); waitbar(i/(Nx-w)); drawnow, end



     % use builtin implementation
     else
         errcode=25;
%         try
         %  X=crp_big(x(i:i+w-1,:),y(i:i+w-1,:),m,t,e,'fan','silent');
           if time_scale_flag
             if length(x(i:i+w-1,:)) > 2000
                Xx=crp_big(x(i:i+w-1,:),m,t,e,method,do_norm{nonorm+1},'silent');
                Xy=crp_big(y(i:i+w-1,:),m,t,e,method,do_norm{nonorm+1},'silent');
                X = double(Xx) .* double(Xy);
             else
                X=jrp(x(i:i+w-1,:),y(i:i+w-1,:),m,t,e,method,do_norm{nonorm+1},'silent');
             end
        else
                Xx=crp2(x(i:i+w-1,:),m,t,e,method,do_norm{nonorm+1},'silent');
                Xy=crp2(y(i:i+w-1,:),m,t,e,method,do_norm{nonorm+1},'silent');
                X = double(Xx) .* double(Xy);
           end
         %  X=crp(x(i:i+w-1,:),y(i:i+w-1,:),m,t,e,varargin{i_char},'silent');

         warning off 
         if nogui~=2 & ishandle(h1), set(h1,'str',[num2str(i),'/',num2str(Nx-w)]); waitbar(i/(Nx-w)); drawnow, end

         if 0
%         catch
           error(lasterr)
           if nogui~=2 & ishandle(hw), close(hw), end
         end

         N=size(X);
         if theiler_window > 0
             X_theiler=double(triu(X,theiler_window))+double(tril(X,-theiler_window));
         else
             X_theiler=X;
         end

         errcode=26;

         % compute recurrence times of 1st and 2nd type
         if size(X_theiler,2) > 1000
             t1=[];t2=[];
             rps2=find(diff(double(X_theiler(:)))==1);
             rps=find(X_theiler(:));
             t1=diff(rps);
             t2=diff(rps2);
         else
             t1 = []; t2 = [];
             for i2=1:size(X_theiler,2)
                 if(Nx-w < 2), waitbar(i2/size(X_theiler,2)), end
                 rps2=find(diff(double(X_theiler(:,i2)))==1);
                 rps=find(X_theiler(:,i2));
                 t1=[t1;diff(rps)];
                 t2=[t2;diff(rps2)];
             end
         end
         t1=mean(t1);
         t2=mean(t2);

         errcode=27;
         [a b]=dl(X_theiler);

         warning off 
         errcode=271;
         b(find(b<lmin))=[];
         [c d]=tt(X_theiler);
         warning off 
         errcode=272;
         d(find(d<vmin))=[];

         errcode=273;
         RR=sum(X_theiler(:))/(N(1)*N(2));
         %b(find(b>=max(N)-lmin))=[]; if isempty(b), b=0; end
         if isempty(b), b=0; end
         errcode=274;
         if sum(sum(X_theiler)) > 0
           DET=sum(b)/sum(sum(X_theiler));
         else
           DET=NaN;
         end
         errcode=275;
         L=mean(b);
         histL=hist(b(:),[1:min(N)]);
         ENTR=entropy(histL(:));
         errcode=276;
         if sum(d)>0
           LAM=sum(d)/sum(sum(X_theiler));
         else
           LAM=NaN;
         end

         errcode=277;
         TT=mean(d);
         b=[b;0]; Lmax=max(b);
         d=[d;0]; Vmax=max(d);
        
    end % end plugin
  
    warning on

    errcode=208;
    Y(i,1)=RR; 
    Y(i,2)=DET;
    Y(i,3)=L;
    Y(i,4)=Lmax;
    Y(i,5)=ENTR;
    Y(i,6)=LAM;
    Y(i,7)=TT;
    Y(i,8)=Vmax;
    Y(i,9)=t1;
    Y(i,10)=t2;
    Y(i,11)=x_var;
    Y(i,12)=y_var;
    Y(i,13)=xy_var;
end
if nogui~=2 & ishandle(hw), waitbar(1); drawnow; close(hw), end

  if ~nogui
    h=findobj('tag','crqa_button_apply');
    set(h(1),'ToolTip','Starts the computation.','String','Apply','Callback','jrqa compute')
    for j=1:length(obj); 
        h=findobj('Tag',obj{j},'Parent',h_fig(1)); 
        if ~isempty(h)
          set(h,'Enable','On')
        end
    end
  end

  if ~nogui, 
    h=findobj('Tag','crqa_Fig'); 
    if ~isempty(h), set(0,'CurrentFigure',h(1)); end
    setptr(gcf,'arrow'), 
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot
if ~nogui
    set(0,'showhidden','on')
    errcode=30;
    h=findobj('Tag','crqa_Fig'); if ~isempty(h), set(0,'CurrentFigure',h(1)); end
    tx={'RR';'DET';'L';'ENTR';'LAM';'TT';'T_1';'T_2';'Variance';'Covariance'};
    index=[1,2,3,5,6,7,9,10,11,12,13];
    tags={'crqa_axes_RR','crqa_axes_DET','crqa_axes_L','crqa_axes_ENTR','crqa_axes_LAM','crqa_axes_TT','crqa_axes_T1','crqa_axes_T2','crqa_axes_Var','crqa_axes_CoVar'};
    h=findobj('Tag','crqa_axes_RR','Parent',gcf); h_axes.h(1)=h(1);
    h=findobj('Tag','crqa_axes_DET','Parent',gcf); h_axes.h(2)=h(1);
    h=findobj('Tag','crqa_axes_L','Parent',gcf); h_axes.h(3)=h(1);
    h=findobj('Tag','crqa_axes_ENTR','Parent',gcf); h_axes.h(4)=h(1);
    h=findobj('Tag','crqa_axes_LAM','Parent',gcf); h_axes.h(5)=h(1);
    h=findobj('Tag','crqa_axes_TT','Parent',gcf); h_axes.h(6)=h(1);
    h=findobj('Tag','crqa_axes_T1','Parent',gcf); h_axes.h(7)=h(1);
    h=findobj('Tag','crqa_axes_T2','Parent',gcf); h_axes.h(8)=h(1);
    h=findobj('Tag','crqa_axes_Var','Parent',gcf); h_axes.h(9)=h(1);
    if ~all(x(:)==y(:)), h=findobj('Tag','crqa_axes_CoVar','Parent',gcf); h_axes.h(10)=h(1); end
   for i=1:9,
     set(gcf,'CurrentAxes',h_axes.h(i))
     if size(Y,1)==1
       cla
       text(0.5,0.5,sprintf('%6.4f',Y(index(i))),'FontWeight','bold','HorizontalAlign','Center')
%        if ~all(x(:)==y(:)) & i==9
%          cla 
%              text(0.5,0.5,sprintf('%6.5f, %6.5f',Y(index(i)),Y(index(i+1))),'FontWeight','bold','HorizontalAlign','Center')
%              set(gcf,'CurrentAxes',h_axes.h(i+1));cla
%          set(h_axes.h(i+1),'visible','off');
%        end
     else
       if i==9
         cla
             h2=stairs(xscale(1:wstep:length(Y)),Y(1:wstep:end,index(i)));
             set(h2,'color',props.line.Color)
         if ~all(x(:)==y(:)) 
           set(gca,'Color','none')
           hold on; h2=stairs(xscale(1:wstep:length(Y)),Y(1:wstep:end,index(10)),'r');
           set(gcf,'CurrentAxes',h_axes.h(i+1));cla
           set(h_axes.h(10),'visible','on');
           h3=stairs(xscale(1:wstep:length(Y)),Y(1:wstep:end,index(11))); set(h3,'color',[0 .4 0]);
           ylabel(tx(i+1));
           set(h_axes.h(i+1),'Tag',tags{i+1},'Units','Norm','Color','none','YAxisLocation','right','YColor',[0 .4 0])
         end
       else
         plot(xscale(1:wstep:length(Y)),Y(1:wstep:end,index(i)),'color',props.line.Color)
       end
     end
     set(gcf,'CurrentAxes',h_axes.h(i));
     ylabel(tx(i));
     set(gca,'Tag',tags{i},'color',props.axes.Color,'Units','Norm')
   end
   h=findobj('Tag','crqa_button_store');
   set(h(1),'UserData',Y,'Enable','On')
   xlim = get(h_axes.h(9),'xlim');
   if length(Y(:,1)) > 1, set(h_crqa_axes_Data,'xlim',xlim), end
   if ~all(x(:)==y(:)) & length(Y(:,1)) == 1, set(h_axes.h(10),'visible','off'); delete(get(h_axes.h(10),'children')), end
else
  if nargout, xout=Y(:,1:10); end
end

if nargout, xout=Y(:,1:10); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the end 

end




if 0
%catch
  if nogui~=2, if ishandle(hw), close(hw), end , end
  
  z=whos;x=lasterr;y=lastwarn;in=varargin{1};
  print_error('jrqa',z,x,y,in,method,action)
  try, if ~nogui
    h=findobj('tag','crqa_button_apply');
    set(h(1),'ToolTip','Starts the computation.','String','Apply','Callback','jrqa compute')
    for j=1:length(obj); 
        h=findobj('Tag',obj{j},'Parent',h_fig(1)); 
        if ~isempty(h)
          set(h,'Enable','On')
        end
    end, end
  end
  set(0,'showhidden','on')
    h=findobj('Tag','crqa_Fig'); if ~isempty(h), setptr(h(1),'arrow'); end
  set(0,'showhidden','off')
end


try, set(0,props.root), end



try
  set(0,'ShowHidden','off')
end
