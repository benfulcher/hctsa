function change_colormapscale(hCRP,cm)
% CHANGE_COLORMAPSCALE   Changes the scale of the current colormap.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 1998-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:36:04 $
% $Revision: 4.8 $
%
% $Log: change_colormapscale.m,v $
% Revision 4.8  2009/03/24 08:36:04  marwan
% copyright address updated
%
% Revision 4.7  2004/11/10 07:04:28  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

  h=findobj('Parent',findobj('Parent',hCRP,'Tag','cm'),'Checked','On');
  c=str2num(get(h,'Tag'));
  h1=cm{c};

  if ~isempty(get(findobj('Tag','CRPData','Parent',findobj('Parent',hCRP,'Tag','CRPPlot')),'UserData'))
    v=[1 2 4 6 8];
%    h1=jet(256);
%    h1=get(hCRP,'Colormap');
    h2=repmat(h1,v(get(findobj('Tag','Log','Parent',gcf),'value')),1);
    h1=h2(1:v(get(findobj('Tag','Log','Parent',gcf),'value')):end,:);
    figure(hCRP)
    colormap(h1)
  end
  clear v h1 h2
