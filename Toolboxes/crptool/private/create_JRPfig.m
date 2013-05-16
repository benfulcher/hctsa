function [h_axes,h_fig]=create_JRPfig(h,xshuttle,yshuttle)
% create_JRPfig   Creates the main figure for the CRP Toolbox
%    Used by CRP Toolbox

% Copyright (c) 2004-2005 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2008/02/25 11:45:27 $
% $Revision: 4.3 $
%
% $Log: create_JRPfig.m,v $
% Revision 4.3  2008/02/25 11:45:27  marwan
% fix of the colorbar bug
%
% Revision 4.2  2005/04/04 09:52:24  marwan
% bug in colormap selection fixed
%
% Revision 4.1  2005/03/16 12:21:30  marwan
% add support for joint recurrence plots
%
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

global props

  scr=get(0,'ScreenSize'); 
  x=xshuttle(:,2:end);
  xscale=xshuttle(:,1);
  y=yshuttle(:,2:end);
  yscale=yshuttle(:,1);

  h_fig=figure('Tag','CRPFig',...		% Plot Figure
            'Position',[scr(3)/4 scr(4)/8 3*scr(3)/8 3*scr(4)/4 ],...
	    'Color',props.window.Color,...
            'NumberTitle','off',...
	    'Name',['Joint Recurrence Plot (' h ')'],...
	    'DeleteFcn','jrp smartclose',...
	    'PaperPosition',[0.25 0.25 7.7677 11.193],...
	    'PaperType','a4',...
	    'PaperOrientation','portrait');
  set(h_fig,props.window,'Units','Norm')
	    
  h1=axes(props.axes,'Units','norm','Position',[.1 .78 .8 .15]);	% Data1 Plot

  if size(y,2)>1
    plot(y,'--','Tag','Data1')
  else
    plot(yscale,y,'r')
  end
  set(h1,'Tag','DataPlot1',...
            'Xcolor','r','ycolor','r',...
	    'XaxisLocation','top',...
	    'YaxisLocation','right',...
	    'UserData',yshuttle)
	    
  h2=axes(props.axes,'Units','norm',...	  % Data2 Plot
           'Position',[.1 .78 .8 .15],...
	    'Xcolor','k','ycolor','k',...
	    'XaxisLocation','bottom','yaxislocation','left');

  if size(y,2)>1
    plot(x,'-','Tag','Data2')
  else
    plot(xscale,x,'k')
  end
  set(h2,'color','none','Tag','DataPlot2','UserData',xshuttle)

  grid on

  if max(abs(x)) > max(abs(y))
      scaling=abs(get(h2,'ylim'));
      if scaling(1)>scaling(2)
         set(h2,'ylim',[-scaling(1) scaling(1)])
      else
         set(h2,'ylim',[-scaling(2) scaling(2)])
      end
      scaling=get(h2,'ylim');
      set(h1,'ylim',scaling)
  else
      scaling=abs(get(h1,'ylim'));
      if scaling(1)>scaling(2)
         set(h1,'ylim',[-scaling(1) scaling(1)])
      else
         set(h1,'ylim',[-scaling(2) scaling(2)])
       end
      scaling=get(h1,'ylim');
      set(h2,'ylim',scaling)
  end

  h1=title('Underlying Time Series','units','normalized');
  h2=get(h1,'Position');
  set(h1,'Position',[h2(1) h2(2)+.12 h2(3)])  

  h_axes=axes(props.axes,'Units','norm',...	% CRP Plot
	    'Color',[1 1 1], ...
            'Tag','CRPPlot',...	
            'Position',[.1 .12 .8 .8*17/23]);
	    
	    
  h1=image( 'Tag','CRPData','cdata',[]);
  h1=title('','units','normalized');
  set(h1,props.titletext)
  h2=get(h1,'Position');
  set(h1,'Position',[h2(1) h2(2)-.03 h2(3)])

  
  h1=text(  .5,.5,'busy...',...			% Text Busy
            props.normaltext,...
	    'Tag','Status',...
            'Visible','off',...
	    'HorizontalAlignment','center',...
	    'VerticalAlignment','middle',...
	    'FontSize',18);

  colormap(french(256))
  h1=colorbar('horiz');				% Colorbar
  set(h1,'Visible','off','Position',[.1 .07 .8 .02],'Tag','Colorbar','xlim',[0 255])
  set(get(h1,'children'),'Visible','off')
  set(findobj('Tag','CRPPlot'),'Position',[.1 .12 .8 .8*17/23])

  cm={'hsv';'hot';'gray';'french';'bone';'copper';...    % Colormap
         'pink';'flag';'lines';'colorcube';...
	 'jet';'prism';'cool';'autumn';...
	 'spring';'winter';'summer'};
  h0=uimenu('Label','Colormap','Tag','cm');
  for i=1:length(cm);
    h1=uimenu(h0,'Label',cm{i},'Checked','Off',...
           'Tag',num2str(i),...
           'Callback','jrp colormap');
    if i==4, set(h1,'Checked','On'), end
  end
  h1=uimenu(h0,'Label','b/w','Tag','18','Callback','jrp colormap');
  h1=uimenu(h0,'Label','inverse','Tag','19','Separator','On','Callback','jrp colormap');

  h1=uimenu('Label','SmartClose',...           % SmartClose
         'Callback','jrp smartclose');

  clear h1 h2 xshuttle yshuttle
  set(h_fig, 'HandleVis','CallBack')
