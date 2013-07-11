function [tout, Nout] = trackplot(varargin)
%TRACKPLOT   Estimates the line of synchronization in a cross recurrence plot. 
%    TRACKPLOT(X [,DX, DY, FP, PARAM]) estimates the line of synchronization in a 
%    cross recurrence plot X. The resulted path can be saved to the 
%    workspace variable T_OUT. This command allows interactive 
%    changing of estimation parameters.
%
%    [A B]=TRACKPLOT(X,DX,DY) estimates the LOS in a recurrence
%    plot X and stores it in A. The number of recurrence points
%    met by the LOS is stored in B(1) and the number of lacks in 
%    the LOS is stored in B(2). 
%
%    TRACKPLOT(X,DX,DY,FP) estimates a LOS whith some fixed points
%    given by the two-column vector FP (eg. [194 0; 201 10]).
%
%    Parameter: The search of the LOS can be forced with the 
%    parameters DX and DY. PARAM can be used to suppress the GUI 
%    (useful in order to use this programme by other programmes).
%
%    Suppressing the GUI.
%      gui         - Creates the GUI and the output plot.
%      nogui       - Suppresses the GUI and the output plot.
%      silent      - Suppresses all output.
%    
%    Examples: y = sin([1:900]*2*pi/67)';
%              y2 = sin(.01*([1:900]*2*pi/67).^2)';
%              x = crp_big(y,y2,3,12,.1,'fan','nogui');
%              trackplot(x,2,2)
%
%    See also CRP2, CRP and CRP_BIG.
% 
%    References:
%    Marwan, N., Thiel, M., Nowaczyk, N.: Cross Recurrence Plot Based 
%    Synchronization of Time Series, Nonlin. Proc. Geophys., 2001.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 1998-2008
% Norbert Marwan, Falko Zetsche, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:33:47 $
% $Revision: 4.12 $
%
% $Log: trackplot.m,v $
% Revision 4.12  2009/03/24 08:33:47  marwan
% copyright address changed
%
% Revision 4.11  2007/12/20 16:26:06  marwan
% changed gpl splash behaviour
%
% Revision 4.10  2005/11/09 08:58:30  marwan
% bug fix in interdependent neighbours method
%
% Revision 4.9  2005/09/02 08:02:57  marwan
% line fitting algorithm improved (linear interpolation between set points)
%
% Revision 4.8  2005/03/16 11:19:02  marwan
% help text modified
%
% Revision 4.7  2004/11/10 07:09:31  marwan
% initial import
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
errcode=0; mflag=0;fixp=[]; X=[];Dmax1=1;Dmax2=1;
set(0,'ShowHidden','On')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check and read the input

warning off
error(nargchk(1,5,nargin))  	% error, if < 1 oder > 4 input variables 
if nargout>2, error('Too many output arguments'), end

if isnumeric(varargin{1})     	% read commandline input
  vargin{4}=[];
  i_double=find(cellfun('isclass',varargin,'double'));
  i_char=find(cellfun('isclass',varargin,'char'));

  % check the text input parameter for gui
  check_gui={'gui','nog','sil'};
  temp_gui=0;
   if ~isempty(i_char)
      for i=1:length(i_char), 
         varargin{i_char(i)}(4)='0';
         temp_gui=temp_gui+strcmpi(varargin{i_char(i)}(1:3),check_gui'); 
      end
      nogui=min(find(temp_gui))-1;
      if isempty(nogui), nogui=0; end
      if nogui>2,nogui=2; end
   else
      nogui=0;
   end
  

   % get the parameters for creating RP
   if max(size(varargin{1}))<2                  
      error('Too less values in data X.')
   end
   X=double(varargin{1});
   Nx = size(X,2); Ny = size(X,1);

   if nargin < 3 ,Dmax1=1; Dmax2=1;
   else
     if isnumeric(varargin{2}), Dmax1=double(varargin{2});
     else Dmax1=1;
     end
     if isnumeric(varargin{3}), Dmax2=double(varargin{3});
     else Dmax2=1;
     end
   end
   if nargin>=4
     if isnumeric(varargin{4}), fixp=fliplr(varargin{4});
     fixp = sortrows(fixp,1);
     else fixp=[];
     end
   end
   
   action='start_gui';

else

%%%%%%%%%%%%%%%%%%%%%%%%%% read input from the GUI

   action=varargin{1};
   nogui=0;
   h=get(gcf,'Name');
   h=h(findstr(h,'(')+1:findstr(h,')')-1);
   hTP=findobj('Name',['TrackPlot (' h ')']);
   hTPCtrl=findobj('Name',['TPControl (' h ')']);
   h=str2num(h);

   temp=get(findobj('Tag','fixp','Parent',findobj('Parent',hTP,'Tag','TPPlot')),'XData');
   if ~isempty(temp)
     for i=1:length(temp), if iscell(temp), fixp(i,1)=temp{i}; else, fixp(i,1)=temp(i); end, end
     temp=get(findobj('Tag','fixp','Parent',findobj('Parent',hTP,'Tag','TPPlot')),'YData');
     for i=1:length(temp), if iscell(temp), fixp(i,2)=temp{i}; else, fixp(i,2)=temp(i); end, end
     fixp = sortrows(fixp,1);
   else 
     fixp=[];
   end
   
   if ( nogui == 0 & ~isempty(hTPCtrl))								
      Dmax1=str2num(get(findobj('Tag','LOSwidthX','Parent',hTPCtrl),'string'));
      Dmax2=str2num(get(findobj('Tag','LOSwidthY','Parent',hTPCtrl),'string'));
      X=double(get(findobj('Tag','RP','Parent',findobj('Parent',hTP,'Tag','TPPlot')),'UserData'));
      Nx=size(X,2); Ny=size(X,1);
   end 


   if isempty(findobj('Tag','TPFig')) & nogui==0 
     action='start_gui';
   end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splash the GPL

splash_gpl('crp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nogui

if nogui>0
   hTP=9999;
   action='LOSsearch';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% switch routines

try
switch(action)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start_gui

case 'start_gui'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% create GUI

  errcode=2;
  if isempty(X)
    lasterr(['??? Error using ==> trackplot',10,...
             'No properly input. Check the class of input.',10,...
	     'Only INT and DOUBLE are allowed, but not LOGICAL.'])
    disp(lasterr)
    return
  end
  scr=get(0,'ScreenSize'); 
  root_ud=get(0,'UserData'); 
  if isstruct(root_ud)
    if isfield(root_ud,'tp')
      if ~isempty(root_ud.tp)
        root_ud.tp=[root_ud.tp max(root_ud.tp)+1];
      else
        root_ud.tp=1;
      end
      h=num2str(root_ud.tp(end));
    else
      root_ud.tp=1;
      h=num2str(1);
    end
  else
    root_ud.old=root_ud;
    root_ud.tp=1;
    h=num2str(1);
  end
  set(0,'UserData',root_ud)

%%%%%%%%%%%%%%%%% Trackplot Figure
   
    h8=figure('Tag','TPFig',...			%Plot Figure
            'Position',[scr(3)/4 scr(4)/4 3*scr(3)/8 3*scr(3)/8 ],...
            'NumberTitle','off',...     
            'Name',['TrackPlot (' h ')'],...
	    'Color',props.window.Color,...
	    'DeleteFcn','trackplot smartclose',...
            'PaperType','a4',...
	    'PaperOrientation','portrait');
    set(h8,props.window,'Units','Norm')

    h1=axes(props.axes,'Parent',h8, ...                    % Initialize CRP Plot
            'Units','normalized', ...
	    'Position',[.1 .1 .8 .8], ...
	    'Color',[1 1 1], ...
	    'Tag','TPPlot', ...
	    'XColor',[0 0 0], ...
	    'YColor',[0 0 0], ...
	    'ZColor',[0 0 0]);

    h2=imagesc('Tag','RP','Parent',h1);               % Plot Data
    X=(X-min(min(X)))/(max(max(X))-min(min(X)));
%    set(h2,'cdata',64*(-(X-max(max(X)))),'UserData',X);
    axis tight
    minvalue=min(min(X));maxvalue=max(max(X));

    set(h1,'tickdir','out','box','on','layer','top')

    h4=line('Parent',h1,'visible','off','Tag','Diagonal',...
       'color',[1 0 0],'LineWidth',1);

  cm={'hsv';'hot';'gray';'french';'bone';'copper';...    % Colormap
         'pink';'flag';'lines';'colorcube';...
	 'jet';'prism';'cool';'autumn';...
	 'spring';'winter';'summer'};
  h0=uimenu('Label','Colormap','Tag','cm');
  if (length(find(X==minvalue))+length(find(X==maxvalue))==length(X(:)))
    set(h2,'cdata',(-X),'UserData',X);
      colormap((gray(2))); cmflag=0;
  else
    set(h2,'cdata',(X),'UserData',X);
      colormap(french(256)); cmflag=1;
  end
  set(h1,'XLim',[0 size(X,2)], 'YLim',[0 size(X,1)])
  for i=1:length(cm);
    h1=uimenu(h0,'Label',cm{i},'Checked','Off',...
           'Tag',num2str(i),...
           'Callback','trackplot colormap');
    if i==18 & cmflag==1, set(h1,'Checked','On'), end
  end
  h1=uimenu(h0,'Label','b/w','Tag','18','Callback','trackplot colormap');
  if cmflag==0, set(h1,'Checked','On'), end
  h1=uimenu(h0,'Label','inverse','Tag','19','Separator','On','Callback','trackplot colormap');

  h1=uimenu('Label','SmartClose',...           % SmartClose
         'Callback','trackplot smartclose');


%%%%%%%%%%%%%%%%% Control Figure

   errcode=3;
   h3=figure('Tag','TPFig',...		% Control Figure
             'Position',[5*scr(3)/8+10 scr(4)/4 1*scr(3)/6 3*scr(4)/9.8],...
             'NumberTitle','off',...
	     'Color',props.window.Color,...
	     'Name',['TPControl (' h ')'],...
	     'DeleteFcn','trackplot handlevisON',...
	     'MenuBar','None',...
	     'Resize','Off');
   set(h3,props.window,'Units','Norm')

   uicontrol(props.frame, ...  % Frame LOSsearch
      	    'Units','Normalized',...
	    'Position',[.1 .51 .8 .42]);


    h0=uicontrol(props.text,...		% Text LOSsearch
            'Units','Normalized',...
	    'FontAngle','italic', ...
	    'String','LOS Search',...
	    'Tag','TextLOSsearch',...
	    'Enable','on',...
	    'Position',[.13 .846 .35 .065]);
    h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h1(3) h1(4)])

    uicontrol(props.button,...		% Set Point
            'Units','Normalized',...
	    'Tag','SetPoint',...
	    'String','Set Point',...
	    'Position',[.16 .765 .314 .08],...
	    'Enable','on',...
	    'ToolTip','Set a Point on LOS.',...
	    'CallBack','trackplot LOSset');

   uicontrol(props.button,...		% Clear Point
            'Units','Normalized',...
	    'Tag','ClearPoint',...
	    'String','Clear P',...
	    'Position',[.52 .765 .314 .08],...
	    'Enable','on',...
	    'ToolTip','Clear a Point on LOS.',...
	    'CallBack','trackplot LOSclear');

   uicontrol(props.button,...		% Clear All Point
            'Units','Normalized',...
	    'Tag','ClearAllPoint',...
	    'String','Clear All',...
	    'Position',[.52 .65 .314 .08],...
	    'Enable','on',...
	    'ToolTip','Clear All Points on LOS.',...
	    'CallBack','trackplot LOSallclear');


   h0=uicontrol(props.text,...			% Text LOSwidthX
            'Units','Normalized',...
	    'Tag','LOSwidthXtext',...
	    'String','dx:',...
	    'Position',[.16 .542 .35 .07],...
	    'Enable','on');
   h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h1(3) h1(4)])

   h0=uicontrol(props.edit,...			% Input LOSwidthX
            'Units','Normalized',...
	    'Tag','LOSwidthX',...
	    'Position',[.31 .55 .14 .06],...
	    'String',num2str(Dmax1),...
	    'Enable','on',...
	    'ToolTip','Insert the LOS search width in X-direction.' );
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h2(3) h1(4)])

  h0=uicontrol(props.text,...			% Text LOSwidthY
            'Units','Normalized',...
	    'Tag','LOSwidthYtext',...
	    'String','dy:',...
	    'Position',[.52 .542 .35 .07],...
	    'Enable','on');
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h1(3) h1(4)])

  h0=uicontrol(props.edit,...			% Input LOSwidthY
            'Units','Normalized',...
	    'Tag','LOSwidthY',...
	    'Position',[.68 .55 .14 .06],...
	    'String',num2str(Dmax2),...
	    'Enable','on',...
	    'ToolTip','Insert the LOS search width in Y-direction.' );
  h1=get(h0,'Extent'); h2=get(h0,'Position'); set(h0,'Position',[h2(1) h2(2) h2(3) h1(4)])


  if ~isunix
    h0=uicontrol(props.frame, ... % Frame Embedding
            'Units','Normalized',...
	    'Position',[.1 .365 .38 .11]);
  end

  uicontrol(props.checkbox, ...       	% Button Stretch Plot
            'Units','Normalized',...
	    'String','Stretch', ...
	    'Position',[.1 .365 .38 .11], ...
	    'Tag','Stretch',...
	    'CallBack','trackplot stretch',...
	    'Value',1,...
	    'ToolTip','Streches the plotted CRP to a squared plot.' );


  uicontrol(props.button,...		% Button Stores LOS
            'Units','Normalized',...
	    'String','Store',...
	    'Position',[.52 .365 .38 .11],...
	    'Tag','Store2',...
	    'Callback','trackplot store2',...
	    'Enable','off',...
	    'ToolTip','Stores the LOS.');

  uicontrol(props.button,...		% Button ApplyLOSsearch
            'Units','Normalized',...
	    'String','Find LOS',...
	    'Position',[.1 .22 .8 .11],...
	    'Tag','Apply2',...
	    'Callback','trackplot LOSsearch',...
	    'Enable','on',...
	    'ToolTip','Searches the LOS.');

  uicontrol(props.button,...		% Button Help
            'Units','Normalized',...
	    'String','Help',...
	    'Position',[.1 .075 .38 .11],...
	    'Tag','Help',...
	    'Callback','helpwin trackplot',...
	    'ToolTip','Opens the helpwindow.');

  uicontrol('Parent',h3, ...
	    props.button,...		% Button Close
            'Units','Normalized',...
	    'String','Close',...
	    'Position',[.52 .075 .38 .11],...
	    'Tag','Close',...
	    'Callback','trackplot close',...
	    'ToolTip','Closes Trackplot windows.');

  set(h8, 'HandleVis','CallBack')
  set(h3, 'HandleVis','CallBack')

  if nargout>=1
    tout=h0;
  end
  if nargout>=2
    Nout=h2;
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change colormap

  case 'colormap'
 
  errcode=81;
  c=str2num(get(gcbo,'Tag'));
  if c~=19,
    set(get(get(gcbo,'Parent'),'Children'),'Checked','Off')
    set(gcbo,'Checked','On')
  end
  cm_old=get(hTP,'Colormap');
  cm=[{hsv(256)}; {hot(256)}; {gray(256)};...
    {french(256)}; {bone(256)}; {copper(256)}; {pink(256)};...
    {flag(256)}; {lines(256)}; {colorcube(256)};...
    {jet(256)};  {prism(256)}; {cool(256)};...
    {autumn(256)}; {spring(256)}; {winter(256)};...
    {summer(256)}; {(gray(2))}; {flipud(cm_old)}];
  set(hTP,'Colormap',(cm{c}))
  
  clear cm cm_old c 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  handlevisON

  case 'handlevisON'
  
  set(hTP, 'HandleVis','on')
	
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% close

  case 'close'  

  errcode=101;
  if ~isempty(findobj('Tag','TPFig')), delete(findobj('Tag','TPFig')), end
  root_ud=get(0,'UserData'); 
  if isstruct(root_ud)
    if isfield(root_ud,'tp')
      root_ud=rmfield(root_ud,'tp');
      if length(fieldnames(root_ud))==1
        if isfield(root_ud,'old'); root_ud=root_ud.old; end
      end
    end
  end
  try, set(0,'UserData',root_ud,props.root), end
  clear all
  disp('Thank you for using CRP toolbox.')


  case 'smartclose'
  errcode=102;

  [h h1]=strtok(get(hTP,'Name'),'(');
  h1([1, end])=[];
  if ishandle(hTPCtrl), delete(hTPCtrl), end
  if ishandle(hTP), delete(hTP), end
  root_ud=get(0,'UserData'); 
  if isstruct(root_ud)
    if isfield(root_ud,'tp')
      root_ud.tp(root_ud.tp==str2num(h1))=[];
      if isempty(root_ud.tp), root_ud=rmfield(root_ud,'tp'); end
    end
    if length(fieldnames(root_ud))==1
      if isfield(root_ud,'old'); root_ud=root_ud.old; end
    end
  end
  try, set(0,'UserData',root_ud,props.root), end
  
  clear all


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% store2

  case 'store2'
  
  errcode=15;
  t=get(findobj('Tag','Apply2','Parent',hTPCtrl),'UserData');
  if isempty(t)
     warndlg('The LOS vector is still empty. Please start the computation of the LOS before storing.','No LOS')
     waitforbuttonpress
  else
     assignin('base','t_out', t)
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stretch

  case 'stretch'

  errcode=7;
  if get(findobj('Tag','Stretch','Parent',gcf),'value')==1
     set(findobj('Tag','TPPlot','Parent',hTP),'PlotBoxAspectRatio',[1 1 1])
  elseif get(findobj('Tag','Stretch','Parent',gcf),'value')==0
     set(findobj('Tag','TPPlot','Parent',hTP),'PlotBoxAspectRatio',[Nx Ny 1])
  end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOSmove

  case 'LOSmove'
  
  errcode=16;
  if isempty(get(gco,'UserData'))
    set(gco,'UserData',1,'ButtonDownFcn','trackplot LOSmove_end')
    set(gcf,'WindowButtonMotionFcn','trackplot LOSmove')
  end
  h1 = round(get(gca,'CurrentPoint'));
  set(gco, 'XData', h1(1,1), 'YData', h1(1,2))
  clear h1

  case 'LOSmove_end'
  if ~isempty(get(gco,'UserData'))
    set(gco,'UserData',[],'ButtonDownFcn','')
    set(gcf,'WindowButtonMotionFcn','')
  end
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOSclear

  case 'LOSclear'
  
  errcode=17;
  if gcf==hTPCtrl
    figure(hTP)
    k=waitforbuttonpress;
    h1 = round(get(gca,'CurrentPoint')); h1(2,:)=[]; h1(:,3)=[];
    finalRect = rbbox;
    h2 = round(get(gca,'CurrentPoint')); h2(2,:)=[]; h2(:,3)=[];
    
    h(1,1)=min(h1(:,1),h2(:,1));h(2,1)=max(h1(:,1),h2(:,1));
    h(1,2)=min(h1(:,2),h2(:,2));h(2,2)=max(h1(:,2),h2(:,2));

    i=find(fixp(:,1)>=h(1,1) & fixp(:,2)>=h(1,2) & fixp(:,1)<=h(2,1) & fixp(:,2)<=h(2,2));
    for j=1:length(i); delete(findobj('tag','fixp','Parent',findobj('Parent',hTP,'Tag','TPPlot'),'xdata',fixp(i(j),1))), end
  else
    h(1,1)=get(gco,'XData');h(1,2)=get(gco,'YData');
    delete(gco)
  end
  clear h

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOSallclear

  case 'LOSallclear'
  
    errcode=171;
    delete(findobj('tag','fixp','Parent',findobj('Parent',hTP,'Tag','TPPlot')))
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOSset

  case 'LOSset'
  
  errcode=18;
  if gcf==hTPCtrl
    figure(hTP)
    ginput(1);
  end
  
  h=round(get(gca,'currentp'));
  fixp(end+1,1)=h(1,1); fixp(end,2)=h(1,2);
  [i j]=sort(fixp(:,1));
  fixp=fixp(j,:);
  h0=uicontextmenu;
  line(h(1,1),h(1,2),1000,...
       'MarkerSize',12,...
       'Marker','.',...
       'Color',[1 0 0],...
       'Tag','fixp',...
       'UIContextMenu',h0)
  uimenu(h0, 'Label', 'Set Point', 'Callback', 'trackplot LOSset')
  uimenu(h0, 'Label', 'Move Point', 'Callback', 'trackplot LOSmove')
  uimenu(h0, 'Label', 'Clear Point', 'Callback', 'trackplot LOSclear')
  clear h h0

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOSsearch

  case 'LOSsearch'
 
  errcode=19;
  x0=0;
  y0=0;
  
  N=size(X);
  Nleer=0; Nvoll=0;

  minvalue=min(min(X));maxvalue=max(max(X));
  if (length(find(X==minvalue))+length(find(X==maxvalue))==length(X(:)))

  flagpoint=0;

% supress output for 'silence'

  if nogui~=2
    h_wait=waitbar(0,'Compute LOS ...'); 
    set(h_wait,'HandleVisibility','on',props.window);
  end
  if nogui==0
    setptr([hTP,hTPCtrl],'watch')
    h=findobj('Tag','Apply2','Parent',hTPCtrl);
    obj.children=get(hTPCtrl,'Children');
    obj.enable=get(obj.children,'Enable'); 
    set(obj.children,'Enable','off')
    set(h(1),'String','Stop',...
             'ToolTip','Stops the computation.',...
             'Enable','on',...
	     'Callback','set(gcbo,''String'',''Stopped'')')
  end

% looks for the beginning of the diagonal LOS

  errcode=191;
  for i=1:N(2);
    if nogui==0
      if check_stop_LOS(hTP,hTPCtrl,nogui,obj), break, end
    end
    if nogui~=2
      waitbar(i/N(2))
    end
    if i<=N(1)
      if ~isempty(fixp), if i>=fixp(1,1), x0=fixp(1,1); y0=fixp(1,2); end, end
      if X(i+y0,1+x0)~=0 
        v=i+y0;
        h=1+x0;
        break
      end
    end
    if X(1+y0,i+x0)~=0
      h=i+x0;
      v=1+y0;
      break
    end
  end
  if ~isempty(fixp) & nogui~=0
    h=fixp(1,1);v=fixp(1,2);x0=0;y0=0;
  end
  Nleer=i-1;
%  t(1:v)=h;
  t(1+x0:h)=v;

% start estimation of the LOS

  errcode=192;
  mflag=0;
  i=2;

  while h<N(2)-1 & v<N(1)-1 & mflag~=1,
   
    if nogui==0
      if check_stop_LOS(hTP,hTPCtrl,nogui,obj), break, end
    end
    if nogui~=2
      waitbar(h/N(2))
    end

    dw=1;
    
    if v < 1, v = 1; end
    if h < 1, h = 1; end
  
    W=X(v:v+dw,h:h+dw);
    W(1,1)=0; 

    if ~isempty(fixp)
      if find(h==fixp(:,1))
        dv=fixp(find(h==fixp(:,1)),2)-v;
        dh=1; 
        if dv <= 0
            h0 = min(find(t > v+dv));
            if isempty(h0), h0 = h; end
            t0 = t(end) + dv;
            fixold = find(h==fixp(:,1));
            if ~isempty(fixold) & fixold > 1; fixold = fixp(fixold - 1,1); else fixold = 0; end
            h1 = max([h0-min([Dmax1,Dmax2]);1;fixold]);
            h2 = max([h0;1;fixold]);
%            h1 = max([h-2*(h-h0);1])
            if length(t) > h1 & length(t) > h2
                t1 = t(h1);
                t2 = t(h2);
                p = polyfit([h1 h2 h],[t1 t2 t0],1);
                t(h1:h) = polyval(p,h1:h);
%                t(h1:h) = spline([h1 h2 h],[t1 t2 t0],h1:h);
            end
            v = v + dv;
            dv = 1; dw = 0;
        end
        flagpoint=1;
      end
    end

    if flagpoint==0;

% looks for the existence of the next recurrence point in diagonal direction

    errcode=193;
    while sum(sum(W))==0 & mflag==0 & flagpoint==0,
      Nleer=Nleer+1;
        if ~isempty(fixp)
           if find(h+dw==fixp(:,1))
  	           W=1; dh=dw;
	           flagpoint=1;
	           break
           end
        end
        dw=dw+1;
        if v+dw < N(1) & h+dw < N(2)
           W=X(v:v+dw,h:h+dw);
           W(1,1)=0; 
        else
           mflag=1;
        end
    end

    if mflag==1
      break
    end

% determines the coordinates of the next recurrence point
  
    errcode=194;
    dh0=min(find(sum(W)));
    dv0=min(find(sum(W')));
    dh=min(find(W(dv0,:)));
    dv=min(find(W(:,dh0)));
    if dh>dh0, dh=dh0;end
    if dv>dv0, dv=dv0;end

% determines the local width of the diagonal LOS

    errcode=195;
    dh1=dh;
    dv1=dv;
      
    % neues Fenster
    WL1=Dmax1;
    WL2=Dmax2;
    if v+dv1-2+WL2>=N(1), WL2=N(1)-(v+dv1-2); end
    if h+dh1-2+WL1>=N(2), WL1=N(2)-(h+dh1-2); end
    Wn=X(v+dv1-1:v+dv1-2+WL2,h+dh-1:h+dh-2+WL1);
    % Schwerpunkt davon ausrechnen
    if sum(sum(Wn))~=0, Sh=sum(sum(Wn).*[1:WL1])/sum(sum(Wn));else Sh=WL1/2;end
    if sum(sum(Wn'))~=0,Sv=sum(sum(Wn').*[1:WL2])/sum(sum(Wn'));else Sv=WL2/2;end
    % neue Dmax berechnen
    if Sh>=Sv, Dmax2n=WL1*Sv/Sh; Dmax1n=WL1; end
    if Sv>Sh,  Dmax1n=WL2*Sh/Sv; Dmax2n=WL2; end

%    while X(v+dv1-1,h:dh-1)==1 & v+dv1-1<N(1) &dv1<Dmax2

    while X(v+dv1-1,h+dh-1)==1 & v+dv1-1<N(1) &dv1<Dmax2n
     dv1=dv1+1;
    end
%    while X(v:dv-1,h+dh1-1)==1 & h+dh1-1<N(2) & dh1<Dmax1
    while X(v+dv-1,h+dh1-1)==1 & h+dh1-1<N(2) & dh1<Dmax1n
      if ~isempty(fixp)
        if find(h+dh1-1==fixp(:,1)) & find(h+fix((dh1+dh)/2)==fixp(:,1))
          dv=fixp(find((h+dh1-1)==fixp(:,1)),2)-v;
	  if isempty(dv), dv=fixp(find((h+fix((dh1+dh)/2))==fixp(:,1)),2)-v; end
          flagpoint=1;
          break
        end
      end
      dh1=dh1+1;
    end


% compute the mean of the diagonal LOS
  
    errcode=196;
    if flagpoint==0
      dh=fix((dh1+dh)/2);
      dv=fix((dv1+dv)/2);
    
      if dh>0
        dh=dh-1;
      end
      if dv>0
        dv=dv-1;
      end
%      dh=dh+dh1;
    end
    
    end % flagpoint end
    flagpoint=0;

% output
    errcode=197;
    if dh~=0 & dv~=0
      t(h:h+dh)=v:dv/dh:v+dv;
    elseif dv~=0
      t(h)=v+dv;
    elseif dv==0
      t(h:h+dh)=v;
    end
    Nvoll=Nvoll+sqrt(dv^2+dh^2);

    if dh==0 & dv==0
      dh=1;
    end

% moves the startpoint for further looking
    h=h+dh;
    v=v+dv;

  end

  else
% DTW algorithm

   if nogui~=2
     h_wait=waitbar(0,'Compute LOS ...');
     set(h_wait,'HandleVisibility','on',props.window);
   end
  
   t=1;
   i=y0+1; j=x0+1;

   while i<N(1)-Dmax1-1 & j<N(2)-Dmax2-1
     errcode=197;
     if nogui~=2
        waitbar(i/N(1)), j0=j;
     end
%     [temp pos]=min([sum(sum(X(i:i+Dmax1,j+1:j+1+Dmax2))),sum(sum(X(i+1:i+1+Dmax1,j+1:j+1+Dmax2))),sum(sum(X(i+1:i+1+Dmax1,j:j+Dmax2)))]);
     [temp pos]=min([(mean(X(i,j+1:j+1+Dmax2))),(mean(diag(X(i+1:i+1+Dmax1,j+1:j+1+Dmax2)))),(mean(X(i+1:i+1+Dmax1,j)))]);
     switch(pos)
       case 1
        j=j+1;
	flag1=1;
       case 2
        i=i+1; j=j+1;
	flag1=1;
	flag2=1;
       case 3
        i=i+1;
	flag2=1;
     end
     if ~isempty(fixp) & flag1==1 & flag2==1
       errcode=1981;
       h=find(fixp(:,1)==j);
       if ~isempty(h) & i<max(fixp(h,2)); 
       h2=find(i<=fixp(h,2)); 
       i=fixp(h(h2(1)),2);  flag1=0; end
     end
     if ~isempty(fixp) & flag2==1 & flag1==1
       errcode=1982;
       h=find(fixp(:,2)==i);
       if ~isempty(h) & j<max(fixp(h,1)); 
       h2=find(j<=fixp(h,1));
       j=fixp(h(h2(1)),1)
       t(j0:j-1)=t(j0);
       flag2=0; end
     end
     t(j)=i;

   end

   t(find(~t))=1;
   Nvoll=NaN; Nleer=NaN;

  end

  errcode=199;
  if nogui~=2
     delete(h_wait)
  end

  if nogui==0
     h1=findobj('Parent',hTP,'Tag','TPPlot');
     h2=findobj('Tag','Diagonal','Parent',h1);
        if isempty(h2)
           h2=line('Parent',h1,'visible','on','Tag','Diagonal','color',[1 0 0],...
          'LineWidth',1);
        end
     set(h2,'xdata',[1:length(t)],'ydata',t,'visible','on')
     h=findobj('Tag','Apply2','Parent',hTPCtrl);
     set(h(1),'String','Apply',...
     	         'ToolTip','Searches the LOS.',...
	         'Callback','trackplot LOSsearch')
     for i=1:length(obj.enable), set(obj.children(i),'Enable',obj.enable{i}); end
     set(findobj('Tag','Store2','Parent',hTPCtrl),'Enable','On') 
     set(findobj('Tag','Apply2','Parent',hTPCtrl),'UserData',t)
  end

  if nargout~=0, tout=t;  end
  if nargout==2
    tout=t; 
    Nout(1)=Nvoll; 
    Nout(2)=Nleer;
  end
  if nogui==0
    setptr([hTP,hTPCtrl],'arrow')
  end

end
warning on
try, set(0,props.root), end
set(0,'ShowHidden','Off')

%%%%%%% error handling

%if 0
catch
  z=whos;x=lasterr;y=lastwarn;in=varargin{1};
  print_error('trackplot',z,x,y,in,[],action)
  try, set(0,props.root), end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=check_stop_LOS(hCRP,hCtrl,nogui,obj)

    global errcode
    errcode=errcode+.02;
    out=0;
    if nogui==0
      h=findobj('Tag','Apply2','Parent',hCtrl);
      if strcmpi(get(h(1),'String'),'stopped')
        set(h(1),'String','Apply',...
     	         'ToolTip','Searches the LOS.',...
	         'Callback','trackplot LOSsearch')
        for i=1:length(obj.enable), set(obj.children(i),'Enable',obj.enable{i}); end
        set(findobj('Tag','Status','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')'),'String','Stopped'),drawnow
        setptr([hCRP,hCtrl],'arrow')
        out=1;
      end
    end
