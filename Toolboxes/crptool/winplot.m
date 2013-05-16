function out=winplot(varargin)
% WINPLOT   Windowed plot
%    WINPLOT(X [,W,WS]) plots means or variances of the
%    sub-vectors of vector X, which have the length W and
%    are shifted by the step WS. X can be a two-column
%    vector, where the first column would be the time-scale.
%
%    WINPLOT(X,W,WS,FLAG) can determine the kind of the
%    result, where FLAG can be either a string or a scalar:
%      'mean' or   1 - Mean (1st moment).
%      'var' or    2 - Variance (2nd moment).
%      'std' or    3 - Standard deviation.
%      'median' or 4 - Median.
%      'sqm' or    5 - Squared Mean.
%      'geo' or    6 - Geometric Mean.
%      '3rd' or    7 - 3rd moment.
%      'skw' or    8 - Skewness.
%      'kur' or    9 - Kurtosis.
%
%    WINPLOT without any arguments calls a demo (the same as the example 
%    below).
%
%    Example: winplot(randn(2000,1),20,20)
%
%    See also PLOT.

% Copyright (c) 2002-2006 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2007/12/20 16:26:07 $
% $Revision: 2.5 $
%
% $Log: winplot.m,v $
% Revision 2.5  2007/12/20 16:26:07  marwan
% changed gpl splash behaviour
%
% Revision 2.4  2006/02/14 11:45:59  marwan
% *** empty log message ***
%
% Revision 2.3  2006/02/07 15:31:41  marwan
% bug in multi-column input resolved
%
% Revision 2.2  2004/11/10 07:04:57  marwan
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check and read the input

slidebar_off=0;

error(nargchk(0,4,nargin));
if nargout>2, error('Too many output arguments.'), end
x=[]; lag=50;

%h=findobj('Tag','errordialog');
%if ~isempty(h), delete(h), end

if isempty(varargin)
  varargin{1}='init';
  x(:,2)=randn(2000,1); 
  x(:,1)=(1:2000)';
  w=20;
  ws=20; 
  flag=1;
end

if ischar(varargin{1})
  action=varargin{1};
else


   i_double=find(cellfun('isclass',varargin,'double'));
   i_char=find(cellfun('isclass',varargin,'char'));

%%%%%% get the numeric input variables
   if ~isempty(i_double)
      for i=1:length(i_double), 
         if i==1
	    x=getx(varargin{i_double(i)});
            w=round(.1*length(x));
            ws=round(.1*length(x));
	    flag=1;
	 elseif i==2
	    w=varargin{i_double(i)};
	    if max(size(w))>1, warning('Window size must be scalar. Using the first value.'); w=w(1); end
	 elseif i==3
	    ws=varargin{i_double(i)};
	    if max(size(ws))>1, warning('Window step value must be scalar. Using the first value.'); ws=ws(1); end
	 elseif i==4
	    flag=varargin{i_double(i)};
	    if max(size(flag))>1, warning('Window step value must be scalar. Using the first value.'); flag=flag(1); end
	    if flag>9, warning(['''',num2str(flag),''' is not a supported method. ''1'' (''mean'') will be used instead.']), flag=1; end
	 end
      end
   else
      error('No data set found.')
   end

%%%%%% get the method
   check_meth={'mea','var','std','med','sqm','geo','3rd','skw','kur'};		% gui, nogui, silent
   if ~isempty(i_char)
      temp_meth=0;
      for i=1:length(i_char), 
         temp_meth=temp_meth+strcmpi(varargin{i_char(i)}(1:3),check_meth'); 
      end
      flag=min(find(temp_meth));
      if isempty(flag) | flag>length(check_meth)
         warning(['''',varargin{i_char(1)},''' is not a supported method. ''mean'' will be used instead.'])
         flag=1;
      end
   end


  if w>length(x), error('Window size cannot be larger than the length of data.'), end

  action='init';

  if nargout
    flag=-flag;
  end

  if flag<0
    try
      out=plot_result(x,w,ws,flag);
      action='none'; 
    catch
      error(lasterr)
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splash the GPL

splash_gpl('crp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% switch routines

try 
switch(action)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initialization
  case 'init'

  if ~isempty(findobj('Tag','winplot_fig'))
%     delete(findobj('Tag','winplot_fig'))
  end


  h=figure(props.window,...				% Plot Figure
            'Tag','winplot_fig',...
	    'MenuBar','Figure',...
            'Position',[69.5000 39.6429 120.0000 30.0714],...
	    'Name','Windowed Plot',...
	    'DeleteFcn','winplot close',...
	    'WindowButtonMotionFcn','winplot motion',...
	    'PaperPosition',[0.25 0.25 7.7677 11.193],...
            'PaperType','a4',...
	    'PaperOrientation','portrait',...
	    'UserData',x);
  
  set(0,'showhidden','on')
  h=findobj('Label','&Help','Type','uimenu');
  if isempty(h)
    h=uimenu('Label','&Help');
    h2=uimenu('Parent',h(1),'Label','&Help Winplot','Callback','helpwin winplot');
  else
    h1=flipud(get(h(1),'Children'));
    set(h1(1),'Separator','on')
    h2=uimenu('Parent',h(1),'Label','&Help Winplot','Callback','helpwin winplot');
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
            'Position',[97 24.2143 22  3.5714]);    
  h2=textwrap(h,{[char(169),' AGNLD'],'University of Potsdam','2002-2006'});
  set(h,'String',h2)

  h=axes(props.axes,...
            'Tag','winplot_axes',...
	    'Box','On',...
            'Position',[9  3.5 72.8333 23.0714]);    
  h=uicontrol(props.frame,...
            'Tag','frame',...
            'Position',[91.5000  1.3571 23.5000 20.7857]);    

 
  h=uicontrol(props.text,...
            'Tag','text',...
	    'String','Window size',...
            'Position',[94.8333 20 16.8333  1.5000]);    

  h=uicontrol(props.edit,...
            'Tag','edit_w',...
	    'String',num2str(w),...
	    'CallBack','winplot update',...
            'Position',[94.8333 18.5714 16.8333  1.5000]);    

  h=uicontrol(props.slider,...
            'Tag','slider_w',...
	    'Max',length(x),...
	    'Min',1,...
	    'Sliderstep',[1/length(x) 10/length(x)],...
	    'Value',w,...
	    'CallBack','winplot wsize',...
            'Position',[94.8333  17.4 16.8333  1.2]);    
  if slidebar_off, set(h,'Visible','off'), end
 
  h=uicontrol(props.text,...
            'Tag','text',...
	    'String','Window Step',...
            'Position',[94.8333 15.2 16.8333  1.5000]);    

  h=uicontrol(props.edit,...
            'Tag','edit_ws',...
	    'String',num2str(ws),...
	    'CallBack','winplot update',...
            'Position',[94.8333 13.7714 16.8333  1.5000]);    

  h=uicontrol(props.slider,...
            'Tag','slider_ws',...
	    'Max',length(x)-1,...
	    'Min',1,...
	    'Value',ws,...
	    'Sliderstep',[1/(length(x)-1) 10/(length(x)-1)],...
	    'CallBack','winplot wstep',...
            'Position',[94.8333  12.4 17.0333  1.2]);    
  if slidebar_off, set(h,'Visible','off'), end

 
  h=uicontrol(props.text,...
            'Tag','text_conf',...
	    'String','Conf.level:',...
            'Position',[94.8333 10.0 16.8333  1.5000]);    

  h=uicontrol(props.edit,...
            'Tag','edit_conf',...
	    'String',num2str(0.05),...
	    'CallBack','winplot update',...
            'Position',[105.6666 10.2 6  1.5000]);    

  h=uicontrol(props.popup,...
            'Tag','button_flag',...
	    'String','Mean|Variance|StdDev|Median|SquaredMean|GeoMean|3rdMoment|Skewness|Kurtosis',...
	    'CallBack','winplot update',...
	    'Value',flag,...
            'Position',[95 7.7857 16.8333  1.8000]);    


  h=uicontrol(props.button,...
            'Tag','button_print',...
	    'CallBack','winplot print',...
	    'String','Print',...
            'Position',[94.8333  4.9286 16.8333  2.2143]);    

  h=uicontrol(props.button,...
            'Tag','button_close',...
	    'CallBack','winplot close',...
	    'String','Close',...
            'Position',[94.8333  2.1429 16.8333  2.2143]);    

  tags={'winplot_fig';'text_conf';'edit_conf';'axes_logo';'text';'winplot_axes';'frame';'edit_w';'slider_w';'edit_ws';'slider_ws';'button_flag';'button_print';'button_close';};
  h=[];
  for i=1:length(tags); h=[h; findobj('Tag',tags{i})]; end
  set(h,'Units','Norm');
  plot_result(x,w,ws,flag);


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wsize
  case 'wsize'

    h=findobj('Tag','slider_w','Parent',gcf);
    w=round(get(h,'Value'));
    h=findobj('Tag','edit_w','Parent',gcf);
    set(h,'String',num2str(w))
    winplot update

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wsstep
  case 'wstep'

    h=findobj('Tag','slider_ws','Parent',gcf);
    ws=round(get(h,'Value'));
    h=findobj('Tag','edit_ws','Parent',gcf);
    set(h,'String',num2str(ws))
    winplot update
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update plot
  case 'update'
    
    h=findobj('Tag','winplot_fig'); h=h(1);
    x=get(h,'UserData');
    h=findobj('Tag','edit_w','Parent',gcf);
    w=round(str2num(get(h,'String')));
    set(h,'String',num2str(w))
    h=findobj('Tag','slider_w','Parent',gcf);
    set(h,'Value',w)
    h=findobj('Tag','edit_ws','Parent',gcf);
    ws=round(str2num(get(h,'String')));
    set(h,'String',num2str(ws))
    h=findobj('Tag','slider_ws','Parent',gcf);
    set(h,'Value',ws)
    h=findobj('Tag','button_flag','Parent',gcf);
    flag=get(h,'Value');
    h=findobj('Tag','edit_conf','Parent',gcf);
    h2=findobj('Tag','text_conf','Parent',gcf);
    if flag~=1 & flag~=2 & flag~=3
      set(h,'Enable','Off')
      set(h2,'Enable','Off')
    else
      set(h,'Enable','On')
      set(h2,'Enable','On')
    end

    plot_result(x,w,ws,flag);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% print
  case 'print'

    h=findobj('Tag','uniLogo');
    h_axes=findobj('Tag','winplot_axes','Parent',gcf); h_axes=h_axes(1);
    h=[h; findobj('Tag','text','Parent',gcf)];
    h=[h; findobj('Tag','frame','Parent',gcf)];
    h=[h; findobj('Tag','edit_w','Parent',gcf)];
    h=[h; findobj('Tag','slider_w','Parent',gcf)];
    h=[h; findobj('Tag','edit_ws','Parent',gcf)];
    h=[h; findobj('Tag','slider_ws','Parent',gcf)];
    h=[h; findobj('Tag','button_flag','Parent',gcf)];
    h=[h; findobj('Tag','button_print','Parent',gcf)];
    h=[h; findobj('Tag','button_close','Parent',gcf)];
    h=[h; findobj('Tag','edit_conf','Parent',gcf)];
    h=[h; findobj('Tag','text_conf','Parent',gcf)];
    set(gcf,'WindowButtonMotionFcn','')
    set(h,'Visible','Off')
    
    set(h_axes,'Units','Character')
    old_pos=get(h_axes,'Position');
    set(h_axes,'Units','normalize','Position',[0.1300 0.1100 0.7750 0.8150])
    h_dlg=printdlg;
    waitfor(h_dlg)

    set(h_axes,'Units','Character','Position',old_pos)
    set(h_axes,'Units','normalize')
    set(h,'Visible','On')
    set(gcf,'WindowButtonMotionFcn','winplot motion')

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% motion
  case 'motion'

    % read data, dim and lag
    h_fig=findobj('Tag','winplot_fig'); h_fig=h_fig(1);
    x=get(h_fig,'UserData');
    h=findobj('Tag','edit_w','Parent',h_fig);
    w=round(str2num(get(h,'String')));
    set(h,'String',num2str(w))
    h=findobj('Tag','slider_w','Parent',h_fig);
    set(h,'Value',w)
    h=findobj('Tag','edit_ws','Parent',h_fig);
    ws=round(str2num(get(h,'String')));
    set(h,'String',num2str(ws))
    h=findobj('Tag','slider_ws','Parent',h_fig);
    set(h,'Value',ws)
    h=findobj('Tag','button_flag','Parent',h_fig);
    flag=get(h,'Value');
    h_axes=findobj('Tag','winplot_axes','Parent',h_fig); 
    h=findobj('Tag','result_line','Parent',h_axes);
    if ~isempty(h)
       y(:,1)=get(h(1),'XData')';
       y(:,2)=get(h(1),'YData')';
    else
       plot_result(x,w,ws,flag);
       h=findobj('Tag','result_line','Parent',h_axes);
       y(:,1)=get(h(1),'XData')';
       y(:,2)=get(h(1),'YData')';
    end
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
    axis_delta=[diff(xlim) diff(ylim)];

    % cursor hover axes?
    pos=get(h_axes,'CurrentPoint');
    if pos(1,1)>=xlim(1) & pos(1,1)<=xlim(2) & pos(1,2)>=ylim(1) & pos(1,2)<=ylim(2)
       % cursor hover line? (within 1% of point centers)
        tolerance=.01;
        is_x1=abs(y(:,1)-pos(1,1)) < tolerance*axis_delta(1);
        is_x2=abs(y(:,1)-pos(2,1)) < tolerance*axis_delta(1);
        is_y1=abs(y(:,2)-pos(1,2)) < tolerance*axis_delta(2);
        is_y2=abs(y(:,2)-pos(2,2)) < tolerance*axis_delta(2);
        res1=is_x1.*is_y1;
        res2=is_x2.*is_y2;

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
	  tx_pos=[y(index1(1),1)+.04*axis_delta(1) y(index1(1),2)+.01*axis_delta(2) 0]; 
	  tx=num2str(y(index1(1),2));
	  set(h,'String',[' ',tx],'Visible','On','Units','Pixel')
	  set(h_index_axes,'Visible','On')
	  tx_ext=get(h,'Extent');
	  set(h_axes,'Units','Pixel')
	  ax_pos=get(h_axes,'Position');
	  set(h_axes,'Units','Norm')
	  xlim=get(h_axes,'Xlim');
	  ylim=get(h_axes,'Ylim');
	  x_offset=(tx_pos(1)-xlim(1))*ax_pos(3)/abs(diff(xlim))+ax_pos(1);
	  y_offset=(tx_pos(2)-ylim(1))*ax_pos(4)/abs(diff(ylim))+ax_pos(2);
	  set(h_index_axes,'Position',[x_offset y_offset tx_ext(3) tx_ext(4)])
          set(h_point,'Position',[y(index1(1),1) y(index1(1),2) 0],'Visible','On')
       elseif ~isempty(index2)
          setptr(h_fig,props.glasspointer)
	  tx_pos=[y(index2(1),1)+.04*axis_delta(1) y(index2(1),2)+.01*axis_delta(2) 0];
	  tx=num2str(y(index2(1),2));
%	  set(h,'Position',tx_pos,'String',[' ',tx],'Visible','On')
          set(h_point,'Position',[y(index2(1),1) y(index2(1),2) 0],'Visible','On')
       else
          set(h_fig,pointer_prefs)
	  set(h,'Visible','Off')
	  set(h_index_axes,'Visible','Off')
       end
    end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% close
  case 'close'

    delete(gcf)

end

try, set(0,props.root), end
set(0,'ShowHidden','off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% error handling
%if 0
catch
  if ~isempty(findobj('Tag','TMWWaitbar')), delete(findobj('Tag','TMWWaitbar')), end
  cmd={['mean'];['var'];['std'];['median'];['squmean'];['geomean'];['bias'];['skewness'];['kurtosis']};
  z=whos;x=lasterr;y=lastwarn;in=varargin{1};
  if ischar(in), in2=in; else, in2=[]; end
  in=whos('in');
  if ~strcmpi(lasterr,'Interrupt')
    fid=fopen('error.log','w');
    err=fprintf(fid,'%s\n','Please send us the following error report. Provide a brief');
    err=fprintf(fid,'%s\n','description of what you were doing when this problem occurred.');
    err=fprintf(fid,'%s\n','E-mail or FAX this information to us at:');
    err=fprintf(fid,'%s\n','    E-mail:  marwan@agnld.uni-potsdam.de');
    err=fprintf(fid,'%s\n','       Fax:  ++49 +331 977 1142');
    err=fprintf(fid,'%s\n\n\n','Thank you for your assistance.');
    err=fprintf(fid,'%s\n',repmat('-',50,1));
    err=fprintf(fid,'%s\n',datestr(now,0));
    err=fprintf(fid,'%s\n',['Matlab ',char(version),' on ',computer]);
    err=fprintf(fid,'%s\n',repmat('-',50,1));
    err=fprintf(fid,'%s\n',x);
    err=fprintf(fid,'%s\n',y);
    err=fprintf(fid,'%s\n',[' during ==> winplot:',action]);
    err=fprintf(fid,'%s\n',[' method ==> ',cmd{flag}]);
    err=fprintf(fid,'%s',[' input ==> ',in.class]);
    if ~isempty(in2), err=fprintf(fid,'\t%s\n',[' (',in2,')']); end
    err=fprintf(fid,'%s\n',[' errorcode ==> no errorcode available']);
    err=fprintf(fid,'%s\n',' workspace dump ==>');
    if ~isempty(z), 
      err=fprintf(fid,'%s\n',['Name',char(9),'Size',char(9),'Bytes',char(9),'Class']);
    for j=1:length(z);
      err=fprintf(fid,'%s\n',[z(j).name,char(9),num2str(z(j).size),char(9),num2str(z(j).bytes),char(9),z(j).class]);
    end, end
    err=fclose(fid);
    disp('------------------------------');
    disp('        ERROR OCCURED');
    disp('   during executing winplot');
    disp('------------------------------');
    disp(x);
    disp(['   during ',action]);
    disp('----------------------------');
    disp('   Please send us the error report. For your convenience, ')
    disp('   this information has been recorded in: ')
    disp(['   ',fullfile(pwd,'error.log')]), disp(' ')
    disp('   Provide a brief description of what you were doing when ')
    disp('   this problem occurred.'), disp(' ')
    disp('   E-mail or FAX this information to us at:')
    disp('       E-mail:  marwan@agnld.uni-potsdam.de')
    disp('          Fax:  ++49 +331 977 1142'), disp(' ')
    disp('   Thank you for your assistance.')
    warning('on')
  end
  set(0,props.root)
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot result

function y=plot_result(x,w,ws,flag)

global props
plotflag=1;
if flag<0; flag=-flag; plotflag=0; end

warning off
if plotflag, setptr(gcf,'watch'), end

cmd={['mean'];['var'];['std'];['median'];['squmean'];['geomean'];['bias'];['skewness'];['kurtosis']};
ylabel_text={'Mean';'Variance';'Std. Deviation';'Median';'Squared Mean';'Geo. Mean';'3rd Moment';'Skewness';'Kurtosis'};

steps=ceil((length(x)+ws-w)/ws);
N=steps*ws+w-ws; x0=x;
if N>length(x), x(end+1:N,:)=NaN; end
if strcmpi(cmd{flag},'geomean')
  if any(any(x<0))
    if plotflag
      h=errordlg('The data have not to be negative for the geometric mean.','Winplot');
      set(h,'Tag','errordialog',props.msgboxwin),h2=guihandles(h);set(h2.OKButton,props.msgbox)
      uiwait(h)
      h=findobj('Tag','button_flag','Parent',gcf);
      set(h,'Value',1)
      setptr(gcf,'arrow')
    else
      error('The data have not to be negative for the geometric mean.')
    end
    return
  end
end

for i=0:steps-1;
  x1=x((i*ws)+1:(i*ws)+w,2);
  x1(isnan(x1))=[]; 
  if isempty(x1), x1=[0;0]; end
  y(i+1)=eval([cmd{flag},'(x1)']);
  vs(i+1)=var(x1);
  vm(i+1)=mean(x1);
  xscale(i+1)=x(i*ws+1,1);
end
xscale=[xscale, x0(end,1)];
y=[y, y(end)];

if plotflag
  h_fig=findobj('Tag','winplot_fig'); h_fig=h_fig(1);
  figure(h_fig(1))
  h_axes=findobj('Tag','winplot_axes','Parent',h_fig); 
  set(h_fig,'CurrentAxes',h_axes)

  h=get(h_axes,'Title');
  old_title=get(h,'String');
  h=get(h_axes,'XLabel');
  old_xlabel=get(h,'String');

  % confidence intervals
  h=findobj('Tag','conf_interval','Parent',h_axes);
  if ~isempty(h)
    delete(h)
  end
  h=findobj('Tag','edit_conf','Parent',gcf);
  alpha=str2num(get(h,'String'));
  n=w-1;

  conf_flag=0;
  switch(flag)

  case 2
    % variance
    conf_up=n*vs/chi2inv(1-alpha/2,n);
    conf_dn=n*vs/chi2inv(alpha/2,n);
    conf_flag=1;

  case 3
    % standard deviation
    conf_up=sqrt(vs/finv(alpha/2,length(x),n));
    conf_dn=sqrt(vs*finv(alpha/2,n,length(x)));
    conf_flag=1;
  
  case 1
    % mean
    conf_up=vm+tinv(1-alpha/2,n)*vs/sqrt(n);
    conf_dn=vm-tinv(1-alpha/2,n)*vs/sqrt(n);
    conf_flag=1;

  end
  if conf_flag
    d_conf_x=diff(xscale); d_conf_x=0;
    conf_x=xscale(1:end-1);
    conf_x=repmat(conf_x,2,1);
    conf_x=conf_x(:); 
    conf_x=[conf_x' 2*conf_x(end)-conf_x(end-1) 2*conf_x(end)-conf_x(end-1) fliplr(conf_x')];
    conf_x(1)=[];
       
    conf_dn=repmat(conf_dn,2,1);
    conf_dn=conf_dn(:); 
    conf_up=repmat(conf_up,2,1);
    conf_up=conf_up(:); 
    conf_y=[conf_dn' fliplr(conf_up') conf_dn(1)];
    
    h=patch(conf_x,conf_y,[1 1 1]);
    set(h,props.patch,'Tag','conf_interval')
    hold on
  end

  % plot stairs line

  h=findobj('Tag','result_line','Parent',h_axes);
  if ~isempty(h)
    delete(h)
  end

  h=stairs(xscale,y);
  set(h,'Tag','result_line',props.line);
  xlabel(old_xlabel)
  title(old_title)

  set(h_axes,'Tag','winplot_axes','YLimMode','auto')

  text(0,0,0,'o','Tag','markedPoint',props.glass)

  grid on
  setptr(gcf,'arrow')
  ylabel(ylabel_text{flag})
  
end
y(2,:)=xscale;
y=rot90(y,-1);
warning on

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% y=getx(x)

function y=getx(x)

if size(x,1)==1
  y(:,1)=(1:length(x))';
  y(:,2)=x(1,:)'; 
else
  if size(x,2)>1 
    if sum(diff(x(:,1))<0)
        y(:,1)=(1:length(x))';
        y(:,2)=x(:,1);
    else
        y(:,1)=x(:,1);
        y(:,2)=x(:,2);
    end
%      error('X-scale must be monotonically non-decreasing.')
  else
    y(:,1)=(1:length(x))';
    y(:,2)=x;
  end
end


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% special function :: squared mean

function y=squmean(x)

y=mean(x).^2;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% special function :: 3rd moment

function y=bias(x)

y=mean((x-mean(x)).^3)^(1/3);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% special function :: skewness

function y=skewness(x)

xs=sort(x); N=length(x);
dz1=xs(ceil(N/10));
dz9=xs(floor(N-N/10));
%q1=xs(ceil(N/4));
%q3=xs(floor(N-N/4));
x2 = 0.33*(dz1 + median(x) + dz9);
y=3*(x2-median(x))/std(x); % Schiefe I
%y=(dz9+dz1-2*median(x))/(dz9-dz1); % Schiefe II
%y=(q3+q1-2*median(x))/(q3-q1); % Schiefe III

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% special function :: kurtosis

function y=kurtosis(x)

% xs=sort(x); N=length(x);
% dz1=xs(ceil(N/10));
% dz9=xs(floor(N-N/10));
% q1=xs(ceil(N/4));
% q3=xs(floor(N-N/4));
% y=(q3-q1)/(2*(dz9-dz1));  % def. Sachs, Angewandte Statistik, 1997
x=(x-mean(x))/std(x);
y=mean(x.^4)-3*mean(x.^2)^2; % def. Hyvaerinen et al., ICA, 2001
