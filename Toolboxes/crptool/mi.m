function varargout=mi(varargin)
%MI   Histogram based mutual information.
%    I=MI(X) computes the auto mutual information using 10 bins,
%    where X can be a multi-column vector.
%
%    I=MI(X1,X2,...,Xn) computes the pair-wise mutual information
%    between all pairs of X1,X2,...,Xn using 10 bins. The result I
%    is a N x N matrix.
%
%    I=MI(X1,X2,...,Xn,N), where N is a scalar, determines the number
%    of bins N.
%
%    I=MI(X1,X2,...,Xn,L) or I=MI(X1,X2,...,Xn,N,L), where L is a 
%    scalar, computes mutual informations for shifting of the 
%    components of each pair of X1,X2,...,Xn until a maximal lag L. 
%
%    [I S]=MI(...) computes the mutual information and the 
%    standard error (only for one- and two-dimensional data).
%
%    MI(...) without any output arguments opens a GUI for interactively 
%    changing the parameters.  
%
%    By using the GUI, the mutual information can be stored into the 
%    workspace. If their standard error is available, they will be
%    appended to the mutual information matrix as the last two columns.
%
%    MI(...,param) additional parameters according to the GUI are
%    available:
%      gui         - Creates the GUI.
%      nogui       - Suppresses the GUI.
%      silent      - Suppresses all output.
%
%    MI without any arguments calls a demo (the same as the example 
%    below).
%
%    Remark
%    Please note that the mutual information derived with MI slightly 
%    differs from the results derived with MIGRAM. The reason is that
%    MI also considers estimation errors. 
%
%    Examples: x = sin(0:.2:8*pi)' + .1*randn(126,1);
%              mi(x,40,100)
%
%    See also HIST2, HISTN, ENTROPY, MIGRAM.
%
%    References: 
%    Roulston, M. S.:
%    Estimating the errors on measured entropy and mutual 
%    information, Physica D, 125, 1999.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2002-2008
% Andre Sitz/ Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2011/05/18 09:46:31 $
% $Revision: 3.16 $
%
% $Log: mi.m,v $
% Revision 3.16  2011/05/18 09:46:31  marwan
% default figure background color bug fixed
%
% Revision 3.15  2010/06/29 12:49:57  marwan
% debugging line removed
%
% Revision 3.14  2009/03/24 08:32:57  marwan
% copyright address changed
%
% Revision 3.13  2007/12/20 16:26:06  marwan
% changed gpl splash behaviour
%
% Revision 3.12  2007/07/17 09:04:38  marwan
% bug for single data series resolved
%
% Revision 3.11  2007/06/19 12:26:28  marwan
% speed enhanced
%
% Revision 3.10  2006/03/29 13:07:55  marwan
% problems regarding OPRPs and embedding resolved
%
% Revision 3.9  2006/02/14 11:46:15  marwan
% *** empty log message ***
%
% Revision 3.8  2005/03/16 11:19:02  marwan
% help text modified
%
% Revision 3.7  2004/11/10 07:05:09  marwan
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check the input

error(nargchk(0,99,nargin));
if nargout>2, error('Too many output arguments'), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splash the GPL

splash_gpl('crp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% error control
%try 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read the input

set(0,'ShowHidden','on')

delete(findobj('Tag','msgbox'))
nogui=0;
lag=[]; 
nbin=10;
x={}; sm_length=[];
if nargin && isnumeric(varargin{1})

  sum_m=0; 

  for i=1:nargin
    if isnumeric(varargin{i}) && max(size(varargin{i}))>1
      y=varargin{i};
      if size(y,1)==1, y=y'; end 
      my=size(y,2); sum_m=sum_m+my;
      mx=size(y,1); 
      if isempty(sm_length), sm_length=mx; end
      if mx < sm_length, sm_length=mx; end
      x(i)={y};
    end
  end
  m=size(x,2); if sum_m==1; lag=0; x(2)=x(1); m=2; end
  if isempty(x), error('Not a valid input vector.'), end

  i_double=find(cellfun('isclass',varargin,'double'));
  i_char=find(cellfun('isclass',varargin,'char'));

  if length(i_double)>1
    if max(size(varargin{i_double(end-1)}))==1
      nbin=varargin{i_double(end-1)}; 
    end
  end
  if max(size(varargin{i_double(end)}))==1
    lag=varargin{i_double(end)}; else if isempty(lag), lag=0; end
  end

  action='init';
  
  if ~isempty(i_char)
    if findstr(lower(varargin{i_char(1)}(1)),'n')
      nogui=1; 
      action='compute'; 
    end
    if findstr(lower(varargin{i_char(1)}(1)),'s')
      nogui=2;
      action='compute'; 
    end
    if findstr(lower(varargin{i_char(1)}(1)),'g')
      nogui=-1;
    end
  end
  if nargout && (nogui<1)
    nogui=1;
    action='compute'; 
  end
  if nogui==-1; nogui=0; end
  	      
elseif nargin && ischar(varargin{1}) && ~isempty(varargin{1})

  if isempty(findobj('Tag','MI_Fig')) && isempty(findobj('Tag','MI_Fig')) 
%    h=errordlg('Start the programme with mi(vector,vector).');
%    set(h,'Tag','msgbox')
    action='end';
  else
  action=varargin{1};
  h=findobj('Tag','MI_Fig');
  x=get(h(1),'UserData');
  h=findobj('Tag','nBin');
  nbin=str2double(get(h(1),'String'));
  m=size(x,2);
  h=findobj('Tag','maxLag');
  lag=str2double(get(h(1),'String'));
  end
else
  action='init';
end

if ~nargin
  x={sin(0:.2:8*pi)'+.1*randn(126,1)}; x(2)=x; m=2;
  nbin=10;
  lag=40; 
  nogui=0;
else
  if ischar(varargin{1})
     action=varargin{1};
  end
end

switch(action)

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% create gui

case 'init'
     
     h8=figure(props.window,...				% Plot Figure
            'Tag','MI_Fig',...
	        'MenuBar','Figure',...
            'Position',[69.5000 36.6429 120.0000 30.0714],...
            'PaperType','a4',...
	        'PaperPosition',[0.25 0.25 7.7677 11.193],...
	        'PaperOrientation','portrait',...
	        'UserData',x,...
  	        'Name','Mutual Information');
  
  set(0,'showhidden','on')
  h=findobj('Label','&Help','Type','uimenu');
  if isempty(h)
    h=uimenu('Label','&Help');
    uimenu('Parent',h(1),'Label','&Help Mutual Information','Callback','helpwin mi');
  else
    h1=flipud(get(h(1),'Children'));
    set(h1(1),'Separator','on')
    uimenu('Parent',h(1),'Label','&Help Mutual Information','Callback','helpwin mi');
    copyobj(h1,h(1))
    delete(h1)
  end

  h=axes(props.axes,...
            'Position',[89 24.8 6.8 3.5]);    
  logo=load('logo');
  h2=imagesc([logo.logo fliplr(logo.logo)]);
  set(h2,'Tag','uniLogo')
  set(h,props.logo,'Tag','axes_logo')
  h=uicontrol(props.text,...
            'Tag','text_logo',...
	        'String','Uni Potsdam',...
            'Position',[97 24.2143 22  3.5714]);    
  h2=textwrap(h,{[char(169),' AGNLD'],'University of Potsdam','2002-2006'});
  set(h,'String',h2)

  axes(props.axes,...
            'Tag','mi_axes',...
	    'Box','On',...
            'Position',[9  3.5 72.8333 23.0714]);    
  uicontrol(props.frame,...
            'Tag','frame',...
            'Position',[91.5000  1.3571 23.5000 20.7857]);    
 
  uicontrol(props.text,...
            'Tag','text',...
	    'String','Number of bins',...
            'Position',[93.8333 20 17.8333  1.5000]);    

  uicontrol(props.edit,...
            'Tag','nBin',...
	    'String',num2str(nbin),...
  	      'ToolTip','Number of bins for estimation of the probability function.',...
            'Position',[93.8333 18.5714 17.8333  1.5000]);    

 
  uicontrol(props.text,...
            'Tag','text',...
	    'String','max. Lag',...
            'Position',[93.8333 16.5 17.8333  1.5000]);    

  uicontrol(props.edit,...
            'Tag','maxLag',...
	    'String',num2str(lag),...
	    'ToolTip','Maximal lag.',...
            'Position',[93.8333 15.0714 17.8333  1.5000]);    

 
   uicontrol(props.text,...
             'Tag','text_show',...
 	     'String','Show components',...
	     'Enable','off',...
             'Position',[93.8333 12.5 21.8333  1.5000]);    
 
   uicontrol(props.popup,...
             'Tag','edit_show1',...
  	     'String','1|2',...
 	     'Enable','off',...
	     'Callback','mi plot_mi',...
             'Position',[93.8333 11.0714  7  1.5000]);    

   uicontrol(props.text,...
             'Tag','text_show',...
 	     'String','vs.',...
	     'Enable','off',...
	     'HorizontalAlign','center',...
             'Position',[101.24995 11.0714-.2 4  1.5000]);    
 
   uicontrol(props.popup,...
             'Tag','edit_show2',...
  	     'String','1|2',...
	     'Callback','mi plot_mi',...
 	     'Enable','off',...
             'Position',[105.6666 11.0714 7 1.5000]);    


   uicontrol(props.button,...
  	      'String','Store',...
  	      'Tag','button_store',...
	      'Enable','Off',...
  	      'ToolTip','Stores the mutual information into a variable in the workspace.',...
  	      'Callback','mi store',...
  	      'Position',[103.8666  8.2143 8.8  2.2143]);

  uicontrol(props.button,...
            'Tag','button_print',...
	    'CallBack','mi print',...
  	      'ToolTip','Prints the MI window.',...
	    'String','Print',...
            'Position',[93.8333  8.2143 8.8  2.2143]);    

  uicontrol(props.button,...
            'Tag','button_close',...
	    'CallBack','mi close',...
  	      'ToolTip','Closes the MI window.',...
	    'String','Close',...
            'Position',[93.8333  2.6429 18.8333  2.2143]);    
  
     uicontrol(props.button,...
  	      'String','Apply',...
  	      'Tag','button_apply',...
  	      'ToolTip','Starts the computation.',...
  	      'Callback','mi compute',...
  	      'Position',[93.8333  5.4286 18.8333  2.2143]);
  
  set(h8, 'HandleVis','CallBack')
  tags={'MI_Fig';'axes_logo';'text_logo';'frame';'text';'mi_axes';'nBin';'maxLag';'text_show';'edit_show1';'edit_show2';'button_store';'button_print';'button_close';'button_apply';};
  h=[];
  for i=1:length(tags); h=[h; findobj('Tag',tags{i})]; end
  set(h,'Units','Norm')
  mi('compute')
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% close windows 

case 'close'
  set(0,props.root)
  h=findobj('Tag','MI_Fig');
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
      MI=get(h(1),'UserData');
      assignin('base',vname, [MI{1:end}])
      h=helpdlg(['Mutual information has been assigned to the workspace variable ''',vname,'''.',10,...
               'If available, the standard errors are given in the second column.'],'Store output');
      set(h,'Tag','msgbox')
      set(h1(1),'UserData',vname)
    end
  end
  set(0,props.root)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% print
  case 'print'

    h=findobj('Tag','uniLogo');
    h_axes=findobj('Tag','mi_axes','Parent',gcf); h_axes=h_axes(1);
    h=[h; findobj('Tag','text','Parent',gcf)];
    h=[h; findobj('Tag','text_logo','Parent',gcf)];
    h=[h; findobj('Tag','frame','Parent',gcf)];
    h=[h; findobj('Tag','nBin','Parent',gcf)];
    h=[h; findobj('Tag','maxLag','Parent',gcf)];
    h=[h; findobj('Tag','text_show','Parent',gcf)];
    h=[h; findobj('Tag','edit_show1','Parent',gcf)];
    h=[h; findobj('Tag','edit_show2','Parent',gcf)];
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

  MI = zeros(m,m,lag);
  MI_sigma = zeros(m,m,lag);

  if ~nogui
    h_fig=findobj('tag','MI_Fig');
    setptr(gcf,'watch'), 
    obj=({'maxLag','nBin','button_store','button_print','button_close','text','text_show','edit_show1','edit_show2'});
    for j=1:length(obj); 
        h=findobj('Tag',obj{j},'Parent',h_fig(1)); 
        if ~isempty(h)
          set(h,'Enable','Off')
        end
    end
    h=findobj('tag','button_apply');
    set(h(1),'ToolTip','Stops the computation.','String','Stop','Callback','set(0,''ShowHidden'',''on'');h=findobj(''tag'',''button_apply'');set(h(1),''String'',''Stopped'');set(0,''ShowHidden'',''off'')')
  end
  if nogui~=2; hw=waitbar(0,'Estimation progress'); end
  for t=0:lag,
     if ~nogui
         set(0,'ShowHidden','on')
         h=findobj('tag','button_apply','Parent',h_fig(1));
         if strcmpi(get(h(1),'string'),'stopped')
             MI((t:lag)+1,1)=NaN;
             MI_sigma((t:lag)+1,1)=NaN;
             break
         end
     end
  
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute the distributions

     for k=1:m; 
         for l=1:m;
             P2=hist2(x{k},x{l},nbin,t,'silent');    % 2D distribution
             mP=length(size(P2));
             %disp(['k: ',num2str(k),' l: ',num2str(l),' test: ',num2str(sum(~P2(:))/numel(P2))])
             if sum(~P2(:))/numel(P2) > 0.9995
                 warning on  
                 warning('Too less data points for the estimation of the joint distribution.'); 
                 MI(k,l,(t:lag)+1)=NaN;
                 MI_sigma(k,l,(t:lag)+1)=NaN;
                 %MI(k,l,(t)+1)=NaN;
                 %MI_sigma(k,l,(t)+1)=NaN;
                 if ~nogui, 
                     h=findobj('tag','button_apply');
                     set(h(1),'ToolTip','Starts the computation.','String','Apply','Callback','mi compute')
                     for j=1:length(obj); 
                         h=findobj('Tag',obj{j},'Parent',h_fig(1)); 
                         if ~isempty(h)
                            set(h,'Enable','On')
                         end
                     end
	                 setptr(h_fig(1),'arrow')
                 end
                 warning on
                 if nogui~=2, delete(hw); end
                 return
             end
             P2=P2/sum(P2(:));      % normalization

             clear P1x
             if mP==2
                 P1x=sum(P2,2);           % distribution for x
                 P1y=sum(P2);             % distribution for y
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute the mutual information
                 I1=[-sum((P1x(P1x~=0)).*log(P1x(P1x~=0))), -sum((P1y(P1y~=0)).*log(P1y(P1y~=0)))];                      % entropies of Px and Py
                 I2=-sum(P2(P2~=0).*log(P2(P2~=0)));                                                                    % entropy of joint distribution Pxy
                 I2_syserr=(length(P1x(P1x~=0))+length(P1y(P1y~=0))-length(P2(P2~=0).*log(P2(P2~=0)))-1)/(2*length(x{k})); % standard error for estimation of entropy
                 MI(k,l,t+1)=I1(1)+I1(2)-I2+I2_syserr;
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute the standard errors
                 P2xy=P1x*P1y;
                 i=(P2xy~=0 & P2~=0);
                 MI_sigma(k,l,t+1)=sqrt( sum( (log(P2xy(i)./P2(i))+MI(t+1)).^2 .*(P2(i).*(1-P2(i)))) /length(x{k}));
             else
                 I1=[];
                 for i=1:mP;
                     P21=permute(P2,[i,i+1:mP, 1:(i-1)]);
                     P1x(:,i)=sum(reshape(P21,size(P21,2),size(P21,2)^(mP-1)),2);   % distribution for x_i
                     I1 = -sum((P1x(P1x~=0)).*log(P1x(P1x~=0))); % sum of entropies of Px_i
                 end  
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute the mutual information
                 I2=-sum(P2(P2~=0).*log(P2(P2~=0)));                                                                    % entropy of joint distribution Pxy
                 I2_syserr=(length(P1x(P1x~=0))-length(P2(P2~=0).*log(P2(P2~=0)))-1)/(2*length(x{k})); % standard error for estimation of entropy
                 MI(k,l,t+1)=I1-I2+I2_syserr;
             end
         end
     end
     
     if nogui~=2; waitbar(t/lag); end
     if ~nogui && t*10/lag == ceil(t*10/lag)
         h=findobj('Tag','button_store');
         if m==2
             out_str{1}=MI; out_str{2}=MI_sigma;
         else
             out_str{1}=MI;
         end
         set(h(1),'UserData',out_str)
         mi('plot_mi')
     end
         
  end
  if nogui~=2; delete(hw); end
  if ~nogui
    set(0,'ShowHidden','on')
    h=findobj('tag','button_apply');
    set(h(1),'ToolTip','Starts the computation.','String','Apply','Callback','mi compute')
    for j=1:length(obj); 
        h=findobj('Tag',obj{j},'Parent',h_fig(1)); 
        if ~isempty(h)
          set(h,'Enable','On')
        end
    end
  end

  if ~isempty(findobj('Tag','MI_Fig'))
    h=findobj('Tag','button_store');
    if m==2
      out_str{1}=MI; out_str{2}=MI_sigma;
    else
      out_str{1}=MI;
    end
    set(h(1),'Enable','On',...
          'UserData',out_str)
  end

  if nargout==1
     varargout(1)={MI};
  elseif nargout==2
     varargout(1)={MI};
     varargout(2)={MI_sigma};
  end
  if ~nogui, setptr(h_fig(1),'arrow'), end
  warning on
  if isnan(MI(1,1)), warning('Mutual information is NaN. Use smaller bin size.'), end
  try set(0,props.root), catch end
  if nargout==0
    mi('plot_mi')
  end
  if nargout==0 && nogui
      varargout{1} = MI;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output

case 'plot_mi'

    h_fig = findobj('Tag','MI_Fig');
    if ~isempty(h_fig)

        h=findobj('Tag','button_store','Parent',h_fig(1));
        if isempty(h), return, end
        x=get(h(1),'UserData');
        MI=x{1}; if m==2; MI_sigma=x{2}; end

        set(0,'current',h_fig(1))
        h_axes=findobj('Tag','mi_axes','Parent',h_fig(1)); h_axes=h_axes(1);
        h_res=findobj('Tag','mi_result','Parent',h_fig(1)); delete(h_res); 
        h_show=[findobj('Tag','text_show','Parent',h_fig(1));...
        findobj('Tag','edit_show1','Parent',h_fig(1));...
        findobj('Tag','edit_show2','Parent',h_fig(1))];
        if lag==0
            set(h_axes,'Visible','off'),cla
	        set(h_show,'enable','off')
            % out_str{MI} = [];
            clear out_str
	        if m==2;
	            for k=1:length(MI)
	                out_str(k)={sprintf(repmat(['%8.3f ',char(177),'%6.3f '],1,2*length(MI)),reshape([MI(k,:);MI_sigma(k,:)],1,2*length(MI)))};
	            end
	        else
	            for k=1:length(MI)
                    out_str(k)={sprintf(repmat('%8.3f ',1,length(MI)),MI(k,:))};
	            end
	        end
            h=uicontrol(props.text,...
                'Tag','mi_result',...
                'HorizontalAlignment', 'left',...
                'FontWeight', 'bold',...
                'String','Mutual Information: ');
            ex=get(h,'Extent'); set(h,'Position',[11  21 ex(3:4)]);
            h=uicontrol(props.listbox,...
                'Tag','mi_result',...
                'Position',[9  7 72.8333 3.0714],...
                'HorizontalAlignment', 'left',...
                'String',out_str);
            ex1=get(h,'Extent'); if ex1(3)>75, ex1(3)=75; end
	        if m<=10
                set(h,'Position',[11  21-m*ex(4) ex1(3) ex1(4)*m]);
	        else
                set(h,'Position',[11  21-10*ex(4) ex1(3) ex1(4)*10]);
	        end

        else
            set(h_axes,'Visible','on'); tx=sprintf('%i|',1:m);
            set(h_show(end-1:end),'String',tx(1:end-1))
            plot(0)
            dimx=get(h_show(end-1),'Value');
            dimy=get(h_show(end),'Value');
            if m==2;
                s1=[permute(MI(dimx,dimy,:),[3,1,2])+3*permute(MI_sigma(dimx,dimy,:),[3,1,2]); flipud(permute(MI(dimx,dimy,:),[3,1,2])-3*permute(MI_sigma(dimx,dimy,:),[3,1,2]))]';
                s3=[0:size(MI,3)-1, size(MI,3)-1:-1:0];
                if all(size(s3) == size(s1))
                    patch(s3,s1,-10,'FaceColor',[.9 .9 1],'EdgeColor',[.85 .85 1])
                end
                hold on
            end
            plot(0:size(MI,3)-1,permute(MI(dimx,dimy,:),[3,1,2])), grid on
            xlabel('Lag'), ylabel('Mutual Information')
            hold off
	        set(gca,'Tag','mi_axes','layer','top','xlim',[0 lag])
        end
        drawnow
        try set(0,props.root), catch end
    end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the end 

case 'end'

end
set(0,'ShowHidden','off')
try set(0,props.root), end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% error handling
if 0
%catch
  if ~isempty(findobj('Tag','TMWWaitbar')), delete(findobj('Tag','TMWWaitbar')), end
  cmd={'mean';'var';'std';'median';'squmean';'geomean';'bias';'skewness';'kurtosis'};
  z_whos=whos;x_lasterr=lasterr;y_lastwarn=lastwarn;if nargin;in=varargin{1};else in.class='no input given';end
  if ischar(in), in2=in; else in2=[]; end
  in=whos('in');
  if ~strcmpi(lasterr,'Interrupt')
    fid=fopen('error.log','w');
    fprintf(fid,'%s\n','Please send us the following error report. Provide a brief');
    fprintf(fid,'%s\n','description of what you were doing when this problem occurred.');
    fprintf(fid,'%s\n','E-mail or FAX this information to us at:');
    fprintf(fid,'%s\n','    E-mail:  marwan@pik-potsdam.de');
    fprintf(fid,'%s\n','       Fax:  ++49 +331 288 2640');
    fprintf(fid,'%s\n\n\n','Thank you for your assistance.');
    fprintf(fid,'%s\n',repmat('-',50,1));
    fprintf(fid,'%s\n',datestr(now,0));
    fprintf(fid,'%s\n',['Matlab ',char(version),' on ',computer]);
    fprintf(fid,'%s\n',repmat('-',50,1));
    fprintf(fid,'%s\n',x_lasterr);
    fprintf(fid,'%s\n',y_lastwarn);
    fprintf(fid,'%s\n',[' during ==> mi:',action]);
    fprintf(fid,'%s',[' input ==> ',in.class]);
    if ~isempty(in2), fprintf(fid,'\t%s\n',[' (',in2,')']); end
    fprintf(fid,'%s\n',' errorcode ==> no errorcode available');
    fprintf(fid,'%s\n',' workspace dump ==>');
    if ~isempty(z_whos), 
        fprintf(fid,'%s\n',['Name',char(9),'Size',char(9),'Bytes',char(9),'Class']);
        for j=1:length(z_whos);
            fprintf(fid,'%s',[z_whos(j).name,char(9),num2str(z_whos(j).size),char(9),num2str(z_whos(j).bytes),char(9),z_whos(j).class]);
            if ~strcmp(z_whos(j).class,'cell') && ~strcmp(z_whos(j).class,'struct')
                  content=eval(z_whos(j).name);
	          content=mat2str(content(1:min([size(content,1),500]),1:min([size(content,2),500])));
                  fprintf(fid,'\t%s',content(1:min([length(content),500])));
            elseif strcmp(z_whos(j).class,'cell')
                  content=eval(z_whos(j).name);
                  fprintf(fid,'\t');
                  for j2=1:min([length(content),500])
                    fprintf(fid,'{%s} ',content{j2});
                  end
            elseif strcmp(z_whos(j).class,'struct')
                  content=fieldnames(eval(z_whos(j).name));
                  content=char(content); content(:,end+1)=' '; content=content';
                  fprintf(fid,'\t%s',content(:)');
            end
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
    disp('----------------------------');
    disp('       ERROR OCCURED');
    disp('    during executing mi');
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
  try 
      if ~nogui
          setptr(h_fig(1),'arrow')
          h=findobj('tag','button_apply');
          set(h(1),'ToolTip','Starts the computation.','String','Apply','Callback','mi compute')
          for j=1:length(obj); 
              h=findobj('Tag',obj{j},'Parent',h_fig(1)); 
              if ~isempty(h)
                set(h,'Enable','On')
              end
          end
      end
  catch
  end
  if nargout, varargout={NaN}; end
  try set(0,props.root),catch end
  set(0,'ShowHidden','off')
end
