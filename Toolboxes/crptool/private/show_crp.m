function show_crp(X,Shuttle)
% SHOW_CRP   Plots the CRP into the GUI
%    Used by CRP Toolbox

% Copyright (c) 1998-2005 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2008/02/25 11:45:27 $
% $Revision: 4.13 $
%
% $Log: show_crp.m,v $
% Revision 4.13  2008/02/25 11:45:27  marwan
% fix of the colorbar bug
%
% Revision 4.12  2007/05/22 13:43:45  marwan
% added fixed RR support
%
% Revision 4.11  2006/03/29 13:07:40  marwan
% problems regarding OPRPs and embedding resolved
%
% Revision 4.10  2005/04/15 09:03:03  marwan
% minor bugfix in plugin section
%
% Revision 4.9  2005/04/08 09:51:46  marwan
% plugin added
%
% Revision 4.8.1.3  2005/04/08 09:03:53  marwan
% plugin added
%
% Revision 4.8.1.2  2005/03/16 12:21:30  marwan
% add support for joint recurrence plots
%
% Revision 4.8.1.1  2004/11/12 08:41:27  marwan
% bug fix in order patterns recurrence plot
%
% Revision 4.8  2004/11/11 12:17:55  marwan
% order patterns recurrence plot added
%
% Revision 4.7  2004/11/10 07:04:29  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

global errcode

NX=size(X,2);
NY=size(X,1);


  switch(Shuttle.mflag)

  case {1,2,3,4,5,6,7,8,9}

  errcode=12;
  set(findobj('Tag','CRPData',...
              'Parent',findobj('Parent',Shuttle.hCRP,'Tag','CRPPlot')),...
              'visible','on',...
              'xdata',Shuttle.xscale(1:NX),...
              'ydata',Shuttle.yscale(1:NY),...
              'cdata',uint8(X),'UserData',X)
  set(Shuttle.hCRP,'colormap',flipud(gray(2)))
  titletext=sprintf('%s \n',[Shuttle.txt_cross,'Recurrence Plot'],...
      ['Dimension: ', num2str(Shuttle.m), ',  Delay: ', num2str(Shuttle.t), ...
       ',  Threshold: ', Shuttle.matext ]);

  if Shuttle.mflag == 8
    titletext=sprintf('%s\n',['Order Matrix']);
  elseif Shuttle.mflag == 9
    titletext=sprintf('%s\n',['Order Patterns ', Shuttle.txt_cross,'Recurrence Plot',10,...
           'Dimension: ', num2str(Shuttle.m), ' (',num2str(factorial(Shuttle.m)),' order patterns),  Delay: ', num2str(Shuttle.t)]);
  end

  case {10,11,12}

  errcode=13;warning on
  X=double(X);
  if ~isempty(get(findobj('Tag','CRPData','Parent',findobj('Parent',Shuttle.hCRP,'Tag','CRPPlot')),'UserData'))
    h1=get(Shuttle.hCRP,'ColorMap');
    if size(h1,1)==2 
       h1=findobj('Parent',findobj('Parent',Shuttle.hCRP,'Tag','cm'),'Checked','On');  
       c=str2num(get(h1,'Tag'));
       v=[1 2 4 6 8];
       h2=repmat(Shuttle.cm{c},v(get(findobj('Tag','Log','Parent',Shuttle.hCtrl),'value')),1);
       h1=h2(1:v(get(findobj('Tag','Log','Parent',Shuttle.hCtrl),'value')):end,:);
    end
    figure(Shuttle.hCRP)
    set(Shuttle.hCRP,'colormap',h1)
  else
    h1=findobj('Parent',findobj('Parent',Shuttle.hCRP,'Tag','cm'),'Checked','On');  
    c=str2num(get(h1,'Tag'));
    h1=Shuttle.cm{c};
    set(Shuttle.hCRP,'colormap',h1)
  end
  zl=0:max(X(:))/9:max(X(:));
  zl=fix(zl*100)/100;  
  X_show=uint8(round(255*X/max(X(:))));
  set(findobj('Tag','CRPData','Parent',findobj('Parent',Shuttle.hCRP,'Tag','CRPPlot')),'visible','on',...
              'xdata',Shuttle.xscale(1:NX),'ydata',Shuttle.yscale(1:NY),...
	      'cdata',X_show,'UserData',X)
  set(get(findobj('Tag','ColBar','Parent',Shuttle.hCRP),'children'),'visible','on')
  
  h1=findobj('Tag','Colorbar','Parent',Shuttle.hCRP);
  set(h1,'visible','on')
  set(get(h1,'xlabel'),'string','Dist. to Next Recurrence Point') 
  b2=get(h1,'position');
  b3=get(h1,'ticklength');

  set(h1,'ticklength',[b3(1)/2 b3(2)],'xlim',[0 255],...
         'xtick',[0:255/9:255],'xticklabel',zl)
  titletext=sprintf('%s \n',[Shuttle.txt_cross,'Distance Matrix'],...
          ['Dimension: ', num2str(Shuttle.m), ',  Delay: ', num2str(Shuttle.t)]);
 
  clear v h1 h2
 
 end % switch
 
  errcode=14;
  set(findobj('Tag','CRPPlot','Parent',Shuttle.hCRP),'tickdir','out','box','on',...
    'xlim',[Shuttle.xscale(1) Shuttle.xscale(NX)],...
    'ylim',[Shuttle.yscale(1) Shuttle.yscale(NY)])

  set(findobj('Tag','DataPlot2','Parent',Shuttle.hCRP),'xlim',[Shuttle.xscale(1) Shuttle.xscale(NX)])
    if get(findobj('Tag','Stretch','Parent',Shuttle.hCtrl),'value')==1
       set(findobj('Tag','CRPPlot','Parent',Shuttle.hCRP),'PlotBoxAspectRatio',[max(NX, NY) max(NX, NY) 1])
       set(findobj('Tag','DataPlot1','Parent',Shuttle.hCRP),'xlim',[Shuttle.yscale(1) Shuttle.yscale(NY)])
    elseif get(findobj('Tag','Stretch','Parent',Shuttle.hCtrl),'value')==0
       set(findobj('Tag','CRPPlot','Parent',Shuttle.hCRP),'PlotBoxAspectRatio',[NX NY 1])
       set(findobj('Tag','DataPlot1','Parent',Shuttle.hCRP),'xlim',[Shuttle.xscale(1) Shuttle.xscale(NX)])
    end

  set(get(findobj('Tag','CRPPlot','Parent',Shuttle.hCRP),'title'),'String',titletext)
  
  set(findobj('Tag','Status','Parent',findobj('Parent',Shuttle.hCRP,'Tag','CRPPlot')),'visible','off',...
         'position',[abs(Shuttle.xscale(1))+abs(Shuttle.xscale(NX)-Shuttle.xscale(1))/2 ...
         abs(Shuttle.yscale(1))+abs(Shuttle.yscale(NY)-Shuttle.yscale(1))/2 0])

  set(findobj('Tag','Status','Parent',findobj('Parent',Shuttle.hCRP,'Tag','CRPPlot')),'String','busy...')
  clear X
