function out = crptool_fnn(varargin)
%FNN   Dimension estimation by false nearest neighbours.
%    Y=FNN(X) computes the vector Y of the amount of false
%    nearest neighbours (FNN) as a function of the embedding dimension.
%
%    Y=FNN(X,M), where M is a scalar, computes the FNN up to dimension M.
%    The defeault is M=10.
%
%    Y=FNN(X,M,T), where T is a scalar, computes the FNN using delay T.
%    The defeault is T=1.
%
%    Y=FNN(X,M,T,R,S), where R and S are scalars, applies the neighbourhood
%    criterion R and the size of the neighbourhood S. The defeault is R=10
%    and S=Inf.
%
%    Y=FNN(X,M,T,R,S,N), where N is a scalar, uses N random samples for
%    the determination of the FNNs. This speeds up the estimation, 
%    especially for long data series. The defeault is N=length(X) if
%    the data length is smaller than 500, else N=200.
%
%    FNN(...) without any output arguments opens a GUI for interactively 
%    changing the parameters. 
%
%    By using the GUI, the FNNs can be stored into the 
%    workspace. 
%
%    FNN(...,param) additional parameters according to the GUI are
%    available:
%      gui         - Creates the GUI.
%      nogui       - Suppresses the GUI.
%      silent      - Suppresses all output.
%
%    FNN without any arguments calls a demo (the same as the example 
%    below).
%
%    Examples: x = sin(0:.2:8*pi)' + .1*randn(126,1);
%              fnn(x,10,[],5)
%
%    See also PHASESPACE, PSS, MI.
%
%    References: 
%    Kennel, M. B., Brown, R., Abarbanel, H. D. I.:
%    Determining embedding dimension for phase-space reconstruction
%    using a geometrical construction, Phys. Rev. A, 45, 1992.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2006-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2014/09/23 07:11:25 $
% $Revision: 5.11 $
%
% $Log: fnn.m,v $
% Revision 5.11  2014/09/23 07:11:25  marwan
% Matlab R2014b incompatibility issues fixed
%
% Revision 5.9  2013/08/06 12:20:41  marwan
% fix of intermixed variables in the code
%
% Revision 5.8  2011/05/18 09:46:17  marwan
% default figure background color bug fixed
%
% Revision 5.7  2010/01/15 12:16:58  marwan
% bug on check of data length
%
% Revision 5.6  2010/01/06 08:59:24  marwan
% bugfix if embedding exceeds data length
%
% Revision 5.5  2009/07/21 11:24:26  marwan
% missing default t
%
% Revision 5.4  2009/03/24 08:35:39  marwan
% copyright address updated
%
% Revision 5.3  2007/12/13 12:02:01  marwan
% added embedding delay
%
% Revision 5.2  2007/05/15 17:33:13  marwan
% new neighbourhood criterion: fixed RR
%
% Revision 5.1  2006/10/24 14:17:47  marwan
% *** empty log message ***
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

global props

init_properties
action='';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check the input

narginchk(0,99);
if nargout>2, error('Too many output arguments'), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splash the GPL

% splash_gpl('crp'); % Do not do this

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% error control
try 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read the input

set(0,'ShowHidden','on')

delete(findobj('Tag','msgbox'))
nogui=0;
maxM_init = 10;
t_init = 1;
maxM=maxM_init; % maximal dimension
r_init = 10;
r=r_init; % r-criterion for neighbours distance
s_init = Inf;
s=s_init; % size of neighbourhood
maxN=[]; % number of samples
x={};
method = [];
if nargin & isnumeric(varargin{1})
    
    % check the numerical input
    i_double=find(cellfun('isclass',varargin,'double'));
    i_length = cellfun('prodofsize',varargin(i_double));
    i_x = i_double(i_length > 1);
    x = varargin{i_x(1)}(:);
    N = length(x);
    t = t_init;
    if isempty(x), error('Not a valid input vector.'), end
    i_x = i_double(i_length < 2);
    switch length(i_x)
      case 1
        maxM = varargin{i_x};
      case 2
        maxM = varargin{i_x(1)};
        t = varargin{i_x(2)};
      case 3
        maxM = varargin{i_x(1)};
        t = varargin{i_x(2)};
        r = varargin{i_x(3)};
      case 4
        maxM = varargin{i_x(1)};
        t = varargin{i_x(2)};
        r = varargin{i_x(3)};
        s = varargin{i_x(4)};
      case 5
        maxM = varargin{i_x(1)};
        t = varargin{i_x(2)};
        r = varargin{i_x(3)};
        s = varargin{i_x(4)};
        maxN = varargin{i_x(5)};
    end

    % check the char input
    i_char=find(cellfun('isclass',varargin,'char'));
    check_meth={'ma','eu','mi'}; 	% maxnorm, euclidean, nrmnorm,  fan, distance
    check_gui={'gui','nog','sil'};				% gui, nogui, silent
    temp_meth=0;
    temp_gui=0;
    
    if ~isempty(i_char)
        for i=1:length(i_char)
            varargin{i_char(i)}(4)='0';
            temp_meth = temp_meth+strcmpi(varargin{i_char(i)}(1:2),check_meth'); 
            temp_gui = temp_gui+strcmpi(varargin{i_char(i)}(1:3),check_gui'); 
        end
    end
    method = min(find(temp_meth));
    nogui = min(find(temp_gui))-1;
    if isempty(nogui); nogui = 0; end


    action = 'init';
    if nogui > 0, action = 'compute'; end
    if nargout & nogui==0
        nogui=1;
        action='compute'; 
    end

elseif nargin & ischar(varargin{1}) & ~isempty(varargin{1})

    if isempty(findobj('Tag','FNN_Fig')) & isempty(findobj('Tag','FNN_Fig')) 
        action='end';
    else

        action=varargin{1};
        h=findobj('Tag','FNN_Fig');
        x=get(h(1),'UserData');
        N=size(x,2);
        h=findobj('Tag','maxM');
        maxM=str2num(get(h(1),'String'));
        h=findobj('Tag','delay');
        t=str2num(get(h(1),'String'));
        h=findobj('Tag','r');
        r=str2num(get(h(1),'String'));
        h=findobj('Tag','s');
        s=str2num(get(h(1),'String'));
        h=findobj('Tag','maxN');
        maxN=str2num(get(h(1),'String'));
        h=findobj('Tag','method');
        method=get(h(1),'Value');
    end
else
    action='init';
end

if ~nargin
      x=sin(0:.2:8*pi)'+.1*randn(126,1); 
      maxM=10; t = 1;
      r=5; s = Inf;
      nogui=0;
else
      if ischar(varargin{1})
          action=varargin{1};
      end
end


% correct invalid inputs
N = length(x);
if isempty(maxM) | maxM > N/2, maxM = maxM_init; end
if isempty(t) | t > N/2 | t < 1, t = t_init; end
if maxM < 2, maxM = 2; end
if (maxM-1) * t > N - 5
    maxM = maxM_init; t = 1; 
    warning('Embedding exceeds data length (m-1) * tau > N. Please check!')
    h = findobj('Tag','delay');
    if ~isempty(h)
        set(h,'string',num2str(t))
    end
    

end
if isempty(r), r = r_init; end
if isempty(s), s = s_init; end
if isempty(maxN)
    maxN = N; 
    if N > 500, maxN = 200; end
end
if isempty(method), method = 1; end

% main part begins
switch(action)

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% create gui

case 'init'
     
    h8=figure(props.window,...				% Plot Figure
            'Tag','FNN_Fig',...
            'MenuBar','Figure',...
            'Position',[69.5000 39.6429 120.0000 30.0714],...
            'PaperType','a4',...
            'PaperPosition',[0.25 0.25 7.7677 11.193],...
            'PaperOrientation','portrait',...
            'UserData',x,...
            'Name','False Nearest Neighbours');
  
    set(0,'showhidden','on')
    h=findobj('Label','&Help','Type','uimenu');
    if isempty(h)
        h=uimenu('Label','&Help');
        h2=uimenu('Parent',h(1),'Label','&Help False Nearest Neighbours','Callback','helpwin fnn');
    else
        h1=flipud(get(h(1),'Children'));
        set(h1(1),'Separator','on')
        h2=uimenu('Parent',h(1),'Label','&Help False Nearest Neighbours','Callback','helpwin fnn');
        copyobj(h1,h(1))
        delete(h1)
    end

    h=axes(props.axes,...
            'Position',[89 24.8+.5 6.8 3.5]);    
    logo=load('logo');
    h2=imagesc([logo.logo fliplr(logo.logo)]);
    set(h2,'Tag','uniLogo')
    set(h,props.logo,'Tag','axes_logo')
    h=uicontrol(props.text,...
            'Tag','text_logo',...
            'String','Uni Potsdam',...
            'Position',[97 24.2143+.5 23  3.5714]);    
    h2=textwrap(h,{[char(169),' AGNLD'],'University of Potsdam','2006'});
    set(h,'String',h2)

    h=axes(props.axes,...
            'Tag','fnn_axes',...
            'Box','On',...
            'Position',[9  3.5 72.8333 23.0714]);    
    
    % frame around parameter settings
    h=uicontrol(props.frame,...
            'Tag','frame',...
            'Position',[91.5000  1.3571 23.5000 22.8]);    

    % dimension
    h=uicontrol(props.text,...
            'Tag','text',...
            'String','Max. Dim.',...
            'Position',[93.8333 22.3-.6 14 1.5000]);    
 
    h=uicontrol(props.edit,...
            'Tag','maxM',...
            'String',num2str(maxM),...
            'ToolTip','False nearest neighbours will be computed up to this dimension.',...
            'Position',[105.1666 22.55-.6 6.5 1.5000]);    

    % delay
    h=uicontrol(props.text,...
            'Tag','text',...
            'String','Delay',...
            'Position',[93.8333 20.3-.6 14 1.5000]);    
 
    h=uicontrol(props.edit,...
            'Tag','delay',...
            'String',num2str(t),...
            'ToolTip','Delay for embedding.',...
            'Position',[105.1666 20.55-.6 6.5 1.5000]);    

    % falseness criterion
    h=uicontrol(props.text,...
            'Tag','text',...
            'String','Falseness',...
            'Position',[93.8333 18.3-.6 20.5  1.5000]);    

    h=uicontrol(props.edit,...
            'Tag','r',...
	        'String',num2str(r),...
	        'ToolTip','Criterion value (suggested: >1).',...
            'Position',[105.1666 18.55-.6 6.5  1.5000]);    
 
    % neigbourhood range
    h=uicontrol(props.text,...
            'Tag','text',...
	        'String','Neigbourh.',...
            'Position',[93.8333 16.3-.6 20.5  1.5000]);    

    h=uicontrol(props.edit,...
            'Tag','s',...
	        'String',num2str(s),...
	        'ToolTip','Maximal distance between neighbours (in multiples of standard deviation).',...
            'Position',[105.1666 16.55-.6 6.5  1.5000]);    
 
    % norm
    h=uicontrol(props.text,...
            'Tag','text',...
	        'String','Norm',...
            'Position',[93.8333 13.7 15  1.5000]);    

    h=uicontrol(props.popup,...
            'Tag','method',...
	        'String','Max|Euc|Min',...
            'Value',method,...
	        'ToolTip','Norm used for the computation of the distance.',...
            'Position',[101.1666 13.95 10.5  1.5000]);    
 
    % number of samples
    h=uicontrol(props.text,...
            'Tag','text',...
	        'String','Number of samples',...
            'Position',[93.8333 11.75+.2 19  1.5000]);    

    h=uicontrol(props.edit,...
            'Tag','maxN',...
	        'String',num2str(maxN),...
	        'ToolTip','Number of randomly chosen samples.',...
            'Position',[93.8333 10.5+.2 17.8333  1.5000]);    
 

    % buttons
    offset = -.5;
    h=uicontrol(props.button,...
        'String','Store',...
        'Tag','button_store',...
        'Enable','Off',...
        'ToolTip','Stores the number of false nearest neighbours into a variable in the workspace.',...
        'Callback','fnn store',...
        'Position',[103.8666  8.2143+offset 8.8  2.2143]);

    h=uicontrol(props.button,...
        'Tag','button_print',...
        'CallBack','fnn print',...
        'ToolTip','Prints the FNN window.',...
        'String','Print',...
        'Position',[93.8333  8.2143+offset 8.8  2.2143]);    

    h=uicontrol(props.button,...
        'Tag','button_close',...
        'CallBack','fnn close',...
        'ToolTip','Closes the FNN window.',...
        'String','Close',...
        'Position',[93.8333  2.6429+offset 18.8333  2.2143]);    

    h=uicontrol(props.button,...
        'String','Apply',...
        'Tag','button_apply',...
        'ToolTip','Starts the computation.',...
        'Callback','fnn compute',...
        'Position',[93.8333  5.4286+offset 18.8333  2.2143]);
  
    set(h8, 'HandleVis','CallBack')
    tags={'FNN_Fig';'axes_logo';'text_logo';'frame';'text';'fnn_axes';'maxM';'r';'s';'method';'maxN';'button_store';'button_print';'button_close';'button_apply';};
    h=[];
    for i=1:length(tags); h=[h; findobj('Tag',tags{i})]; end
    set(h,'Units','Norm')
    fnn('compute')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% close windows 

case 'close'
  set(0,props.root)
  h=findobj('Tag','FNN_Fig');
  if ~isempty(h), close(h(1)), end 
  clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% store

case 'store'
  if ~isempty(findobj('Tag','button_store'))
    h=findobj('Tag','button_store');
    h1=findobj('Tag','button_close');
    if ~isempty(h1), vname_old=get(h1(1),'UserData'); else vname_old=''; end
    if isempty(vname_old), vname_old=''; end
    vname=char(inputdlg('Choose a variable name.','Store output',1,{vname_old}));
    if isempty(vname)
      return
    else
      FNN=get(h(1),'UserData');
      assignin('base',vname, [FNN])
      h=helpdlg(['False nearest neighbours have been assigned to the workspace variable ''',vname,'''.',...
               ''],'Store output');
      set(h,'Tag','msgbox')
      set(h1(1),'UserData',vname)
    end
  end
  set(0,props.root)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% print
  case 'print'

    h=findobj('Tag','uniLogo');
    h_axes=findobj('Tag','fnn_axes','Parent',gcf); h_axes=h_axes(1);
    h=[h; findobj('Tag','text','Parent',gcf)];
    h=[h; findobj('Tag','text_logo','Parent',gcf)];
    h=[h; findobj('Tag','frame','Parent',gcf)];
    h=[h; findobj('Tag','maxM','Parent',gcf)];
    h=[h; findobj('Tag','r','Parent',gcf)];
    h=[h; findobj('Tag','s','Parent',gcf)];
    h=[h; findobj('Tag','method','Parent',gcf)];
    h=[h; findobj('Tag','maxN','Parent',gcf)];
    h=[h; findobj('Tag','button_store','Parent',gcf)];
    h=[h; findobj('Tag','button_print','Parent',gcf)];
    h=[h; findobj('Tag','button_close','Parent',gcf)];
    h=[h; findobj('Tag','button_apply','Parent',gcf)];
    set(h,'Visible','Off')
    
    set(h_axes,'Units','Character')
    old_pos=get(h_axes,'Position');
    set(h_axes,'Units','normalize','Position',[0.1300 0.1100 0.7750 0.8150])
    h_dlg=printdlg;
    waitfor(h_dlg)

    set(h_axes,'Units','Character','Position',old_pos)
    set(h_axes,'Units','normalize')
    set(h,'Visible','On')
    set(0,props.root)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% computation 

case 'compute'
    warning off
    if ~nogui
        h_fig=findobj('tag','FNN_Fig');
        set(0,'ShowHidden','on')
        h_axes=findobj('Tag','fnn_axes','Parent',gcf); h_axes=h_axes(1);
        h_stopbutton = findobj('tag','button_apply','Parent',h_fig(1));
        h_storebutton = findobj('Tag','button_store');
        setptr(gcf,'watch'), 
        obj=({'maxM','method','r','s','maxN','button_store','button_print','button_close','text'});
        for j=1:length(obj); 
            h=findobj('Tag',obj{j},'Parent',h_fig(1)); 
            if ~isempty(h), set(h,'Enable','Off'), end
        end
        set(h_stopbutton(1),'ToolTip','Stops the computation.','String','Stop','Callback','set(0,''ShowHidden'',''on'');h=findobj(''tag'',''button_apply'');set(h(1),''String'',''Stopped'');set(0,''ShowHidden'',''off'')')
    end
    if nogui~=2; 
        hw=waitbar(0,'Estimation progress'); 
%        hw=waitbar(0,'Estimation progress','CreateCancelBtn','return'); 
    end
   

    FNN(1,1) = 1;
    x_s = s * std(x);
    for m = 2:maxM
        % check if stop button was pressed
        if ~nogui
            set(0,'ShowHidden','on')
            if strcmpi(get(h_stopbutton(1),'string'),'stopped')
                FNN((m:maxM)+1,1)=NaN;
                break
            end
        end
        
        FNN(m,1) = 0;
        NX=N-t*(m-1);

        % create phase space vectors
        clear jx,for mi=1:m;
            jx(1+NX*(mi-1):NX+NX*(mi-1))=1+t*(mi-1):NX+t*(mi-1);
        end

        x1=reshape(repmat(x,1,m),N*m,1);
        x2=reshape(x1(jx),NX,m); % phasespace vector
        jx = reshape(jx,NX,m);
        if maxN >= NX; 
            idx = 1:NX;
        else
            idx = ceil(NX * rand(maxN,1)); %random samples
        end
     
        % count the false neighbours
        cnt = 0;
        for i = 1:length(idx), if nogui~=2 & i*10/length(idx) == ceil(i*10/length(idx)); waitbar((length(idx)*(m-1)+i)/(length(idx)*maxM)), end
            % distance in the m-dimensional space
            switch method
              case 1
                distance = max(abs(x2-repmat(x(jx(idx(i),:))',size(x2,1),1))');
              case 2
                distance = sqrt(sum((x2-repmat(x(jx(idx(i),:))',size(x2,1),1)).^2,2));
              case 3
                distance = sum(abs(x2-repmat(x(jx(idx(i),:))',size(x2,1),1)),2);
            end
            [distance is] = sort(distance);
            
            if (is(2)+m*t) < N & (idx(i)+m*t) < N & distance(2) < x_s
                %FNN(m,1) = FNN(m,1) + double((abs(x(idx(i)+m) - x(is(2)+m))/distance(2) > r));
                FNN(m,1) = FNN(m,1) + double((abs(x(idx(i)+(m)*t) - x(is(2)+(m)*t))/distance(2) > r));
                cnt = cnt + 1;
            end
            % check if stop button was pressed
            if ~nogui
                if strcmpi(get(h_stopbutton(1),'string'),'stopped')
                    FNN((m:maxM)+1,1)=NaN;
                    break
                end
            end
        end
        FNN(m,1) = FNN(m,1) / cnt;
        
        if ~nogui, 
            %fnn('plot_fnn')
            plot(h_axes,1:length(FNN),FNN)
            xlabel(h_axes,'Dimension'), ylabel(h_axes,'False Nearest Neighbours')
            set(h_axes,'Tag','fnn_axes','layer','top','xlim',[1 maxM],'ylim',[0 1],'xtick',[1:1:maxM],'xgrid','on','ygrid','on')
            drawnow
        end
    
    end

    if nogui~=2 & ishandle(hw), delete(hw); end
    if ~nogui, 
        set(0,'ShowHidden','on')
        h=findobj('tag','button_apply');
        set(h(1),'ToolTip','Starts the computation.','String','Apply','Callback','fnn compute')
        for j=1:length(obj); 
            h=findobj('Tag',obj{j},'Parent',h_fig(1)); 
            if ~isempty(h), set(h,'Enable','On'), end
        end
	    setptr(h_fig(1),'arrow')
        if ~isempty(findobj('Tag','FNN_Fig'))
                out_str=FNN;
                set(h_storebutton(1),...
                      'UserData',out_str)
            h=findobj('Tag','button_store');
            set(h(1),'Enable','On')
        end

        if ~nogui, setptr(h_fig(1),'arrow'), end
        warning on
        try, set(0,props.root), end
        if nargout==0
            fnn('plot_fnn')
        end
    end
    if nogui
       out=FNN;
    end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output

case 'plot_fnn'
    set(0,'ShowHidden','on')
    h_fig=findobj('tag','FNN_Fig');
    if ~isempty(h_fig)
        h=findobj('Tag','button_store','Parent',h_fig(1));
        if isempty(h), return, end
        FNN=get(h(1),'UserData');

        h_axes=findobj('Tag','fnn_axes','Parent',h_fig(1)); h_axes=h_axes(1);
        h_res=findobj('Tag','fnn_result','Parent',h_fig(1)); delete(h_res); 
        set(0,'current',h_fig(1))
 
        set(h_axes,'Visible','on');
        plot(1:length(FNN),FNN), grid on
        xlabel('Dimension'), ylabel('False Nearest Neighbours')
        set(gca,'Tag','fnn_axes','layer','top','xlim',[1 maxM],'ylim',[0 1],'xtick',[1:1:maxM])
        drawnow
    end
    try, set(0,props.root), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the end 

case 'end'

end
set(0,'ShowHidden','off')
try set(0,props.root), end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% error handling
%if 0
catch
  if ~isempty(findobj('Tag','TMWWaitbar')), delete(findobj('Tag','TMWWaitbar')), end
  cmd={['mean'];['var'];['std'];['median'];['squmean'];['geomean'];['bias'];['skewness'];['kurtosis']};
  z_whos=whos;x_lasterr=lasterr;y_lastwarn=lastwarn;if nargin;in=varargin{1};else in.class='no input given';end
  if ischar(in), in2=in; else, in2=[]; end
  in=whos('in');
  if ~strcmpi(lasterr,'Interrupt')
    fid=fopen('error.log','w');
    err=fprintf(fid,'%s\n','Please send us the following error report. Provide a brief');
    err=fprintf(fid,'%s\n','description of what you were doing when this problem occurred.');
    err=fprintf(fid,'%s\n','E-mail or FAX this information to us at:');
    err=fprintf(fid,'%s\n','    E-mail:  marwan@pik-potsdam.de');
    err=fprintf(fid,'%s\n','       Fax:  ++49 +331 288 2640');
    err=fprintf(fid,'%s\n\n\n','Thank you for your assistance.');
    err=fprintf(fid,'%s\n',repmat('-',50,1));
    err=fprintf(fid,'%s\n',datestr(now,0));
    err=fprintf(fid,'%s\n',['Matlab ',char(version),' on ',computer]);
    err=fprintf(fid,'%s\n',repmat('-',50,1));
    err=fprintf(fid,'%s\n',x_lasterr);
    err=fprintf(fid,'%s\n',y_lastwarn);
    err=fprintf(fid,'%s\n',[' during ==> fnn:',action]);
    err=fprintf(fid,'%s',[' input ==> ',in.class]);
    if ~isempty(in2), err=fprintf(fid,'\t%s\n',[' (',in2,')']); end
    err=fprintf(fid,'%s\n',[' errorcode ==> no errorcode available']);
    err=fprintf(fid,'%s\n',' workspace dump ==>');
    if ~isempty(z_whos), 
      err=fprintf(fid,'%s\n',['Name',char(9),'Size',char(9),'Bytes',char(9),'Class']);
    for j=1:length(z_whos);
      err=fprintf(fid,'%s',[z_whos(j).name,char(9),num2str(z_whos(j).size),char(9),num2str(z_whos(j).bytes),char(9),z_whos(j).class]);
      if ~strcmp(z_whos(j).class,'cell') & ~strcmp(z_whos(j).class,'struct')
            content=eval(z_whos(j).name);
	    content=mat2str(content(1:min([size(content,1),500]),1:min([size(content,2),500])));
            err=fprintf(fid,'\t%s',content(1:min([length(content),500])));
      elseif strcmp(z_whos(j).class,'cell')
            content=eval(z_whos(j).name);
            err=fprintf(fid,'\t');
            for j2=1:min([length(content),500])
	      content2=content{j2};
              err=fprintf(fid,'{%s} ',content{j2});
            end
      elseif strcmp(z_whos(j).class,'struct')
            content=fieldnames(eval(z_whos(j).name));
            content=char(content); content(:,end+1)=' '; content=content';
            err=fprintf(fid,'\t%s',content(:)');
      end
      err=fprintf(fid,'\n');
    end, end
    err=fclose(fid);
    disp('----------------------------');
    disp('       ERROR OCCURED');
    disp('    during executing fnn');
    disp('----------------------------');
    disp(x_lasterr);
    disp(['   during ',action]);
    disp('----------------------------');
    disp('   Please send us the error report. For your convenience, ')
    disp('   this information has been recorded in: ')
    disp(['   ',fullfile(pwd,'error.log')]), disp(' ')
    disp('   Provide a brief description of what you were doing when ')
    disp('   this problem occurred.'), disp(' ')
    disp('   E-mail or FAX this information to us at:')
    disp('       E-mail:  marwan@pik-potsdam.de')
    disp('          Fax:  ++49 +331 288 2640'), disp(' ')
    disp('   Thank you for your assistance.')
    warning('on')
  end
  try if ~nogui
    setptr(h_fig(1),'arrow')
    h=findobj('tag','button_apply');
    set(h(1),'ToolTip','Starts the computation.','String','Apply','Callback','fnn compute')
    for j=1:length(obj)
        h=findobj('Tag',obj{j},'Parent',h_fig(1)); 
        if ~isempty(h)
          set(h,'Enable','On')
        end
    end
  end,end
  if nargout, out=NaN; end
  try set(0,props.root), end
  set(0,'ShowHidden','off')
end
