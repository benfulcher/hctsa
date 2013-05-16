function phasespace(varargin)
% PHASESPACE   Embedding of time series in a phase space
%    PHASESPACE(X) shows the 3D phase space trajectory of 
%    the system which is presented by the observation X.
%    The phase vectors are a reconstruction by using the
%    time delay method (Takens, 1981). A GUI provides to
%    change the embedding dimension to 2D.
%
%    PHASESPACE(X,Y) uses the two one-column vectors
%    X and Y as components of the phase space trajectory.
%    The representation is 2D only and cannot be switched
%    to 3D.
%
%    PHASESPACE(X,Y,Z) uses the three one-column vectors
%    X, Y and Z as components of the phase space trajectory.
%    The representation is 3D only and cannot be switched
%    to 2D.
%
%    PHASESPACE without any arguments calls a demo (the
%    same as the example below).
%
%    Example: phasespace(cos(0:.1:32).*[321:-1:1])
%
%    See also FNN, PSS.
%
%    References: 
%    Takens, F.: 
%    Detecting Strange Attractors in Turbulence,
%    Lecture Notes in Mathematics, 898, Springer, Berlin, 1981

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2002-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:34:24 $
% $Revision: 2.5 $
%
% $Log: phasespace.m,v $
% Revision 2.5  2009/03/24 08:34:24  marwan
% copyright address changed
%
% Revision 2.4  2007/12/20 16:26:06  marwan
% changed gpl splash behaviour
%
% Revision 2.3  2006/03/29 13:07:55  marwan
% problems regarding OPRPs and embedding resolved
%
% Revision 2.2  2006/02/14 11:46:15  marwan
% *** empty log message ***
%
% Revision 2.1  2004/11/10 07:07:56  marwan
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

global props

init_properties
try

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check and read the input

slidebar_off=0;

error(nargchk(0,3,nargin));
x=[]; lag=50;

if isempty(varargin)
  varargin{1}='init';
end
if ischar(varargin{1})

action=varargin{1};

else

if nargin==3
  dim=3;
  x(:,3)=varargin{3}(:,1);
  x(:,2)=varargin{2}(:,1);
  x(:,1)=varargin{1}(:,1);
elseif nargin==2
  dim=2;
  x(:,2)=varargin{2}(:,1);
  x(:,1)=varargin{1}(:,1);
elseif nargin==1
  dim=3;
  x=varargin{1};
  if size(x,1)==1, x=x(1,:)'; else, x=x(:,1); end
end

m=size(x,2);
if m>3
  h=errordlg('Only the first three columns of the vector will be used.','Vector size');
  waitfor(h);
  m=3; x(:,4:end)=[];
end


action='init';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splash the GPL

splash_gpl('crp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% switch routines

switch(action)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initialization
  case 'init'

  if isempty(x), x(:,1)=(cos(0:.1:62).*[621:-1:1])'; m=1; dim=3; lag=53; end

  h=figure(props.window,...				% Plot Figure
            'Tag','phasespace_fig',...
	    'MenuBar','Figure',...
            'Position',[69.5000 39.6429 120.0000 30.0714],...
	    'Name','Phase Space Reconstruction',...
	    'DeleteFcn','phasespace close',...
	    'WindowButtonMotionFcn','phasespace motion',...
	    'PaperPosition',[0.25 0.25 7.7677 11.193],...
            'PaperType','a4',...
	    'PaperOrientation','portrait',...
	    'UserData',x);
  ha=guihandles(h);
  ha_fn=fieldnames(ha);
  for i=1:length(ha);
    try
      set(eval(['ha.',ha_fn{i}]),props.menubar)
    end
  end
  
  set(0,'showhidden','on')
  h=findobj('Label','&Help','Type','uimenu');
  if isempty(h)
    h=uimenu('Label','&Help');
    h2=uimenu('Parent',h(1),'Label','&Help Phasespace','Callback','helpwin phasespace');
  else
    h1=flipud(get(h(1),'Children'));
    set(h1(1),'Separator','on')
    h2=uimenu('Parent',h(1),'Label','&Help Phasespace','Callback','helpwin phasespace');
    copyobj(h1,h(1))
    delete(h1)
  end
  set(0,'showhidden','off')

  h=axes(props.axes,...
            'Position',[89 24.8 6.8 3.5]);    
  logo=load('logo');
  h2=imagesc([logo.logo fliplr(logo.logo)]);
  set(h2,'Tag','uniLogo')
  set(h,props.logo,'Tag','axes_logo')
  h=uicontrol(props.text,...
            'Tag','text',...
	    'String','Uni Potsdam',...
  	    'Units','Character',...
            'Position',[97 24.2143 22  3.5714]);    
  h2=textwrap(h,{[char(169),' AGNLD'],'University of Potsdam','2002-2006'});
  set(h,'String',h2)

  h=axes(props.axes,...
            'Tag','phasespace_axes',...
	    'Box','On',...
            'Position',[4.8333  3.5000 76.8333 25.0714]);    
  h=uicontrol(props.frame,...
            'Tag','frame',...
            'Position',[91.5000  1.3571 23.5000 20.7857]);    


  h=uicontrol(props.slider,...
            'Tag','slider_hori',...
	    'CallBack','phasespace rotate',...
	    'Min',-360,...
	    'Max',360,...
	    'SliderStep',[1/360 10/360],...
            'Position',[4.8333  0.6429 76.8333  1.3500]);    
  h=uicontrol(props.slider,...
            'Tag','slider_vert',...
	    'CallBack','phasespace rotate',...
	    'Min',-90,...
	    'Max',90,...
	    'SliderStep',[1/180 10/180],...
            'Position',[85  3.5000  3.000 25.0714]);    


  h=uicontrol(props.text,...
	    'String','Dimension',...
            'Tag','text',...
            'Position',[94.8333 20.1429 16.8333  1.5000]);    
  if m~=1; set(h,'Enable','Off'), end	    
  h=uicontrol(props.radio,...
            'Tag','button_2d',...
	    'String','2D',...
	    'Interruptible','off',...
	    'CallBack','phasespace dim2',...
	    'Value',0,...
            'Position',[94.8333 18.8857 16.8333  1.5000]);    
  if dim==2, set(h,'Value',1), end
  if m~=1; set(h,'Enable','Off'), end	    
  h=uicontrol(props.radio',...
            'Tag','button_3d',...
	    'String','3D',...
 	    'Value',0,...
	    'Interruptible','off',...
	    'CallBack','phasespace dim3',...
            'Position',[94.8333 17.3857 16.8333  1.5000]);    
  if dim==3, set(h,'Value',1), end
  if m~=1; set(h,'Enable','Off'), end	    
 
  h=uicontrol(props.text,...
            'Tag','text',...
	    'String','Lag',...
            'Position',[94.8333 15.4 16.8333  1.5000]);    
  if m~=1; set(h,'Enable','Off'), end	    
  if lag>(.1*length(x)), lag=ceil(.1*length(x)); end
  h=uicontrol(props.edit,...
            'Tag','edit_lag',...
	    'String',num2str(lag),...
	    'CallBack','phasespace update',...
	    'Interruptible','off',...
            'Position',[94.8333 14.1214 16.8333  1.5000]);    
  if m~=1; set(h,'Enable','Off'), end	    
  h=uicontrol(props.slider,...
            'Tag','slider_lag',...
	    'Max',length(x),...
	    'Interruptible','off',...
	    'Value',lag,...
	    'Min',0,...
	    'Sliderstep',[1/length(x) 10/length(x)],...
	    'CallBack','phasespace lag',...
            'Position',[94.8333  12.95 16.8333  1.2]);    
  if m~=1; set(h,'Enable','Off'), end	    
  if slidebar_off, set(h,'Visible','off'), end

  h=uicontrol(props.radio,...
            'Tag','button_line',...
	    'Interruptible','off',...
	    'String','Lines',...
	    'CallBack','phasespace linestyle',...
	    'Value',1,...
            'Position',[94.8333 10.7857 16.8333  1.5000]);

  h=uicontrol(props.radio,...
            'Tag','button_dots',...
	    'String','Dots',...
	    'CallBack','phasespace linestyle',...
	    'Interruptible','off',...
	    'Value',0,...
            'Position',[94.8333 9.2857 16.8333  1.5000]);

  h=uicontrol(props.radio,...
            'Tag','button_comet',...
	    'String','Comet',...
	    'Interruptible','off',...
	    'CallBack','phasespace linestyle',...
	    'Value',0,...
            'Position',[94.8333 7.7857 16.8333  1.5000]);


  h=uicontrol(props.button,...
            'Tag','button_print',...
	    'CallBack','phasespace print',...
	    'String','Print',...
            'Position',[94.8333  4.9286 16.8333  2.2143]);    

  h=uicontrol(props.button,...
            'Tag','button_close',...
	    'CallBack','phasespace close',...
	    'String','Close',...
            'Position',[94.8333  2.1429 16.8333  2.2143]);    

  tags={'phasespace_fig';'axes_logo';'text';'slider_hori';'slider_vert';'phasespace_axes';'frame';'button_2d';'button_3d';'edit_lag';'slider_lag';'button_line';'button_dots';'button_comet';'button_print';'button_close';};
  h=[];
  for i=1:length(tags); h=[h; findobj('Tag',tags{i})]; end
  set(h,'Units','Norm');
  plot_trajectory(x,dim,lag)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dimension
  case 'dim2'

    h=findobj('Tag','button_2d','Style','RadioButton','Parent',gcf);
    if get(h,'Value')~=1
      h=findobj('Tag','button_2d','Parent',gcf); 
      set(h(1),'Value',1);
    end
    h=findobj('Tag','button_3d','Parent',gcf);
    set(h(1),'Value',0);
    phasespace update

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dimension
  case 'dim3'

    h=findobj('Tag','button_3d','Style','RadioButton','Parent',gcf); 
    if get(h,'Value')~=1
      h=findobj('Tag','button_3d','Parent',gcf);
      set(h(1),'Value',1);
    end
    h=findobj('Tag','button_2d','Parent',gcf);
    set(h(1),'Value',0);
    phasespace update

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% lag (slider)
  case 'lag'

    h=findobj('Tag','slider_lag','Style','Slider','Parent',gcf); 
    lag=round(get(h,'Value'));
    h=findobj('Tag','edit_lag','Parent',gcf);
    set(h(1),'String',num2str(lag))
    phasespace update
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% linestyle
  case 'linestyle'

    h_line=findobj('Tag','button_line','Parent',gcf); h_line=h_line(1);
    h_dots=findobj('Tag','button_dots','Parent',gcf); h_dots=h_dots(1);
    h_comet=findobj('Tag','button_comet','Parent',gcf); h_comet=h_comet(1);
    h=findobj('Tag','phasespace_trajectory'); h=h(1);
    flag=get(h_comet,'Value');
    if h_line==gcbo
      set(h,props.line)
      set(h_line,'Value',1)
      set(h_dots,'Value',0)
      set(h_comet,'Value',0)
    elseif h_dots==gcbo
      set(h,props.marker)
      set(h_line,'Value',0)
      set(h_dots,'Value',1)
      set(h_comet,'Value',0)
    elseif h_comet==gcbo
      set(h_line,'Value',0)
      set(h_dots,'Value',0)
      set(h_comet,'Value',1)
    end
    if flag
      phasespace update
    end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update trajectory
  case 'update'
    
    h=findobj('Tag','phasespace_fig'); h=h(1);
    x=get(h,'UserData');
    h=findobj('Tag','edit_lag','Parent',gcf); 
    lag=round(str2num(get(h,'String')));
    set(h,'String',num2str(lag))
    h=findobj('Tag','slider_lag','Parent',gcf);
    set(h,'Value',lag)
    h=findobj('Tag','button_3d','Parent',gcf); 
    if get(h,'Value')
      dim=3;
    end
    h=findobj('Tag','button_2d','Parent',gcf); 
    if get(h,'Value')
      dim=2;
    end

    plot_trajectory(x,dim,lag)
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% rotate axes
  case 'rotate'

    h=findobj('Tag','slider_vert','Parent',gcf); 
    el=get(h,'value');
    h=findobj('Tag','slider_hori','Parent',gcf); 
    az=get(h,'value');
    view(az,el)
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% print
  case 'print'

    h=findobj('Tag','uniLogo');
    h_axes=findobj('Tag','phasespace_axes','Parent',gcf); h_axes=h_axes(1);
    h=[h; findobj('Tag','text','Parent',gcf)];
    h=[h; findobj('Tag','frame','Parent',gcf)];
    h=[h; findobj('Tag','slider_hori','Parent',gcf)];
    h=[h; findobj('Tag','slider_vert','Parent',gcf)];
    h=[h; findobj('Tag','button_2d','Parent',gcf)];
    h=[h; findobj('Tag','button_3d','Parent',gcf)];
    h=[h; findobj('Tag','edit_lag','Parent',gcf)];
    h=[h; findobj('Tag','slider_lag','Parent',gcf)];
    h=[h; findobj('Tag','button_line','Parent',gcf)];
    h=[h; findobj('Tag','button_dots','Parent',gcf)];
    h=[h; findobj('Tag','button_comet','Parent',gcf)];
    h=[h; findobj('Tag','button_print','Parent',gcf)];
    h=[h; findobj('Tag','button_close','Parent',gcf)];
    set(gcf,'WindowButtonMotionFcn','')
    set(h,'Visible','Off')
    
    old_pos=get(h_axes,'Position');
    set(h_axes,'Units','normalize','Position',[0.1300 0.1100 0.7750 0.8150])
    h_dlg=printdlg;
    waitfor(h_dlg)

    set(h_axes,'Position',old_pos)
    set(h,'Visible','On')
    set(gcf,'WindowButtonMotionFcn','phasespace motion')

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% motion
  case 'motion'

  h_fig=findobj('Tag','phasespace_fig'); h_fig=h_fig(1);

  h=findobj('Tag','button_comet','Style','RadioButton','Parent',h_fig);
  if get(h,'Value')~=1
    % read data, dim and lag
    h=findobj('Tag','phasespace_fig'); h=h(1);
    x=get(h,'UserData');
    m=size(x,2);
    h=findobj('Tag','edit_lag','Style','Edit','Parent',h_fig);
    lag=round(str2num(get(h,'String')));
    set(h,'String',num2str(lag))
    h=findobj('Tag','slider_lag');
    set(h(1),'Value',lag)
    h=findobj('Tag','button_3d','Style','RadioButton','Parent',h_fig);
    if get(h,'Value')
      dim=3;
    end
    h=findobj('Tag','button_2d','Style','RadioButton','Parent',h_fig);
    if get(h,'Value')
      dim=2;
    end
    h2=findobj('Tag','phasespace_trajectory','Type','Line');
    if isempty(h2)
       plot_trajectory(x,dim,lag);
    end

    % create vectors
    switch(m)
      case {2,3}
        y=x;
	if dim==2, y(:,3)=0; end
      case 1
        if dim==3, y=[x(1:end-2*lag),x(1+lag:end-lag),x(1+2*lag:end)]; 
        elseif dim==2, y=[x(1:end-lag),x(1+lag:end)]; y(:,3)=0;
        end
    end

    h_axes=findobj('Tag','phasespace_axes','Parent',h_fig); h_axes=h_axes(1);
    pointer_prefs=getptr(h_fig);
    pointer_prefs=cell2struct(pointer_prefs(2:2:end),pointer_prefs(1:2:end),2);
    bad_pointer=setptr(props.glasspointer);
    bad_pointer=cell2struct(bad_pointer(2:2:end),bad_pointer(1:2:end),2);
    fnames=fieldnames(pointer_prefs);
    i=find(strcmpi(fnames,'Pointer'));
    if strcmpi(getfield(pointer_prefs,fnames{i}),'custom')
      i=find(strcmpi(fnames,'PointerShapeCData'));
      shape_pointer=getfield(pointer_prefs,fnames{i});
      shape_pointer(isnan(shape_pointer))=0;
      fnames_bad_pointer=fieldnames(bad_pointer);
      i=find(strcmpi(fnames_bad_pointer,'PointerShapeCData'));
      shape_bad_pointer=getfield(bad_pointer,fnames_bad_pointer{i});
      shape_bad_pointer(isnan(shape_bad_pointer))=0;
      if ~prod(prod(double(shape_pointer==shape_bad_pointer)))
        set(h_axes,'UserData',pointer_prefs)
      else
        pointer_prefs=get(h_axes,'UserData');
      end
    else
      set(h_axes,'UserData',pointer_prefs)
    end

    xlim=get(h_axes,'xlim');
    ylim=get(h_axes,'ylim');
    zlim=get(h_axes,'zlim');
    axis_delta=[diff(xlim) diff(ylim) diff(zlim)];

    % cursor hover axes?
    pos=get(h_axes,'CurrentPoint');
    if dim==2, pos(:,3)=0; end
    if (pos(1,1)>=xlim(1) & pos(1,1)<=xlim(2) & pos(1,2)>=ylim(1) & pos(1,2)<=ylim(2)  & pos(1,3)>=zlim(1) & pos(1,3)<=zlim(2)) | ...
       (pos(2,1)>=xlim(1) & pos(2,1)<=xlim(2) & pos(2,2)>=ylim(1) & pos(2,2)<=ylim(2)  & pos(2,3)>=zlim(1) & pos(1,3)<=zlim(2))
       % cursor hover line? (within 1% of point centers)
        if dim==2, tolerance=.01; else tolerance=.05; end
        is_x1=abs(y(:,1)-pos(1,1)) < tolerance*axis_delta(1);
        is_x2=abs(y(:,1)-pos(2,1)) < tolerance*axis_delta(1);
        is_y1=abs(y(:,2)-pos(1,2)) < tolerance*axis_delta(2);
        is_y2=abs(y(:,2)-pos(2,2)) < tolerance*axis_delta(2);
%         is_x1=(y(:,1)>pos(1,1) & y(:,1)<pos(2,1)) + (y(:,1)<pos(1,1) & y(:,1)>pos(2,1));
%         is_y1=(y(:,2)>pos(1,2) & y(:,2)<pos(2,2)) + (y(:,2)<pos(1,2) & y(:,2)>pos(2,2));
        res1=is_x1.*is_y1;
        res2=is_x2.*is_y2;
%          is_z1=(y(:,3)>pos(1,3) & y(:,3)<pos(2,3)) + (y(:,3)<pos(1,3) & y(:,3)>pos(2,3));
        is_z1=abs(y(:,3)-pos(1,3)) < tolerance*axis_delta(3);
        is_z2=abs(y(:,3)-pos(2,3)) < tolerance*axis_delta(3);
        res1=res1.*is_z1;
        res2=res2.*is_z2;

       % show index and change cursor
       index1=find(res1);index2=find(res2);
       h_index_axes=findobj('Tag','index_axes','Parent',h_fig);
       h=findobj('Tag','index_text','Parent',h_index_axes);
       if isempty(h_index_axes)
          h_index_axes=axes(props.axesindex,'Tag','index_axes');
       end  
       if isempty(h)
          h=text(props.textindex,'Units','norm','Tag','index_text','Parent',h_index_axes);
       end  
       h=h(1); h_index_axes=h_index_axes(1);
       h_point=findobj('Tag','markedPoint','Parent',h_axes); 
       set(h_point,'Visible','Off')
       if ~isempty(index1)
          setptr(h_fig,props.glasspointer)
	  tx_pos=[y(index1(1),1)+.04*axis_delta(1) y(index1(1),2)+.01*axis_delta(2) y(index1(1),3)+.01*axis_delta(3)]; 
	  tx=num2str(index1(1));
	  set(h,'String',[' ',tx],'Visible','On','Units','Pixel')
	  set(h_index_axes,'Visible','On')
	  tx_ext=get(h,'Extent');
	  set(h_fig,'Units','Pixel')
          pos=get(h_fig,'CurrentPoint');
	  set(h_fig,'Units','Norm')
	  set(h_index_axes,'Position',[pos(1)+props.index.xoffset pos(2)+props.index.yoffset tx_ext(3) tx_ext(4)])
          set(h_point,'Position',[y(index1(1),1) y(index1(1),2) y(index1(1),3)],'Visible','On')
       elseif ~isempty(index2)
          setptr(h_fig,props.glasspointer)
	  tx_pos=[y(index2(1),1)+.04*axis_delta(1) y(index2(1),2)+.01*axis_delta(2) y(index2(1),3)+.01*axis_delta(3)];
	  tx=num2str(index2(1));
	  set(h,'String',[' ',tx],'Visible','On','Units','Pixel')
	  set(h_index_axes,'Visible','On')
	  tx_ext=get(h,'Extent');
	  set(h_fig,'Units','Pixel')
          pos=get(h_fig,'CurrentPoint');
	  set(h_fig,'Units','Norm')
	  set(h_index_axes,'Position',[pos(1)+props.index.xoffset pos(2)+props.index.yoffset tx_ext(3) tx_ext(4)])
          set(h_point,'Position',[y(index2(1),1) y(index2(1),2) y(index2(1),3)],'Visible','On')
       else
          set(h_fig,pointer_prefs)
	  set(h,'Visible','Off')
 	  set(h_index_axes,'Visible','Off')
      end
    end
  end   

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% close
  case 'close'

   delete(gcf)
   
end

try, set(0,props.root), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% error handling
catch
%if 0
  z=whos;x=lasterr;y=lastwarn;in=varargin{1};
  errcode=[];
  print_error('phasespace',z,x,y,in,[],action)
end



 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot trajectory

function plot_trajectory(x,dim,lag)

global props

try
if (length(x)-(dim-1)*lag)<=0
  h=warndlg(['No phasespace vectors available.',char(10),...
            'Please use a smaller lag.'],'Warning'); 
  set(h,props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
  set(h,'Tag','helpdlgbox','HandleVisibility','On');
  return
end

h_fig=findobj('Tag','phasespace_fig'); h_fig=h_fig(1);
figure(h_fig)
h_axes=findobj('Tag','phasespace_axes','Parent',h_fig); h_axes=h_axes(1);
h=findobj('Tag','button_comet','Parent',h_fig); h=h(1);
if get(h,'Value')==1
  cmd='comet';
else
  cmd='plot';
end
flag='';

m=size(x,2);
switch(m)
  case 3
    parm='x(:,1),x(:,2),x(:,3)';
    flag='3';
  case 2
    parm='x(:,1),x(:,2)';
  case 1
    if dim==3
      parm='x(1:end-2*lag),x(1+lag:end-lag),x(1+2*lag:end)';
      flag='3';
    elseif dim==2
      parm='x(1:end-lag),x(1+lag:end)';
    end
end

set(h_fig,'CurrentAxes',h_axes)
eval([cmd,flag,'(',parm,')']);
set(h_axes,'Tag','phasespace_axes')
h=findobj('Type','Line','Parent',gca);
h_line=findobj('Tag','button_line','Parent',h_fig); h_line=h_line(1);
h_dots=findobj('Tag','button_dots','Parent',h_fig); h_dots=h_dots(1);
if get(h_line,'Value')==1
  set(h,props.line)
elseif get(h_dots,'Value')==1
  set(h,props.marker)
end
set(h,'Tag','phasespace_trajectory')

if strcmp(cmd,'plot')
hold on
% marker points in trajectory
switch(m)
  case 3
    parm='x(i,1),x(i,2),x(i,3)';
  case 2
    parm='x(i,1),x(i,2)';
  case 1
    if dim==3, parm='x(i),x(i+lag),x(i+2*lag)'; 
    elseif dim==2, parm='x(i),x(i+lag)';
    end
end

m_size=[30, 25, 20, 15, 10];
dx=abs(round(.001*diff(get(h_axes,'xlim'))/mean(diff(x(:,1)))));
for j=1:5, 
  i0=round(.25*length(x-(dim-1)*lag));
  i=i0+dx*j;
  if i<(length(x)-(dim-1)*lag)
    eval(['plot',flag,'(',parm,',''.'',''markersize'',',num2str(m_size(j)),')'])
  end
  i0=round(.75*length(x-(dim-1)*lag));
  i=i0+dx*j;
  if i<(length(x)-(dim-1)*lag)
    eval(['plot',flag,'(',parm,',''.'',''markersize'',',num2str(m_size(j)),')'])
  end
end
text(0,0,0,'o','Tag','markedPoint',props.glass)
hold off
axis equal
end

[az, el]=view;
h=findobj('Tag','slider_hori','Parent',gcf);
set(h,'Value',az)
h=findobj('Tag','slider_vert','Parent',gcf);
set(h,'Value',el)
if strcmp(cmd,'plot')
  grid on
end
catch 
end
