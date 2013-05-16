function x_out=choosecolormap(x_in)
%CHOOSECOLORMAP   GUI for choosing a colormap
%   CHOOSECOLORMAP enables to change the colormap of the 
%   current figure.
%
%   See also COLORMAP, GRAPH3D

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2002-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:30:23 $
% $Revision: 2.2 $
%
% $Log: choosecolormap.m,v $
% Revision 2.2  2009/03/24 08:30:23  marwan
% copyright address changed
%
% Revision 2.1  2004/11/10 07:07:20  marwan
% initial import
%
%

error(nargchk(0,1,nargin)) 
if nargin==0
  x_in=gcf;
end

% defining the colormap
cm.name={'hsv';'hot';'gray';'bone';'copper';...
         'pink';'flag';'lines';'colorcube';...
	 'vga';'jet';'prism';'cool';'autumn';...
	 'spring';'winter';'summer'};
cm.value(1)={hsv};
cm.value(2)={hot};
cm.value(3)={gray};
cm.value(4)={bone};
cm.value(5)={copper};
cm.value(6)={pink};
cm.value(7)={flag};
cm.value(8)={lines};
cm.value(9)={colorcube};
cm.value(10)={vga};
cm.value(11)={jet};
cm.value(12)={prism};
cm.value(13)={cool};
cm.value(14)={autumn};
cm.value(15)={spring};
cm.value(16)={winter};
cm.value(17)={summer};


% make the GUI
if ~ischar(x_in)
  h=get(0,'Children');
  if max(h==x_in)==0
    error('Specified figure does not exist.')
  end
  
  h0=figure('Tag','Choosecolormap',...
             'NumberTitle','off',...
	     'Name',['Colormap'],...
	     'MenuBar','None',...
	     'UserData',x_in,...
	     'Units','Char');
  hp=get(h0,'Position'); set(h0,'Position',[hp(1) hp(2) 30 20])
  h1=uicontrol('Style','Listbox',...
               'Tag','CMchoice',...
               'Units','Norm',...
               'Position',[0 0 1 1],...
	       'String',cm.name,...
	       'Callback','choosecolormap apply');

else
  switch(x_in)
  case 'apply'
    fig=get(gcf,'UserData');
    choice=get(findobj('Tag','CMchoice'),'Value');
    set(fig,'ColorMap',cm.value{choice})
    close gcf
  end
end
