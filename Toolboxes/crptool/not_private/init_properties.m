function init_properties
% INIT_PROPERTIES   Initializes the basic properties for GUIs
%    Used by CRP Toolbox

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
% $Log: init_properties.m,v $
% Revision 4.8  2009/03/24 08:36:04  marwan
% copyright address updated
%
% Revision 4.7  2004/11/10 07:04:29  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

global props

    font.FontName=get(0,'factoryTextFontName');
    font.FontUnits=get(0,'factoryTextFontUnits');
    font.FontSize=get(0,'factoryTextFontSize');
    font.FontWeight=get(0,'factoryTextFontWeight');
    font.FontAngle=get(0,'factoryTextFontAngle');
    props.root.defaultUIControlBackgroundColor=get(0,'factoryUIControlBackgroundColor');
    props.root.defaultUIControlForegroundColor=get(0,'factoryUIControlForegroundColor');
    props.root.defaultTextColor=get(0,'factoryUIControlForegroundColor');
    props.window.Color=[0.8008 0.7500 0.6875];
    props.window.Units='Char';
    props.window.Colormap=[[0 0 0];[props.window.Color]];
    props.menubar.Color=props.window.Color;
    props.frame=font;
    props.frame.Style='frame';    
    props.frame.BackgroundColor=props.window.Color;
    props.frame.ForegroundColor=[0 0 0];
    props.frame.Units='Char';
    props.slider=font;
    props.slider.Style='slider';    
    props.slider.BackgroundColor=props.window.Color;
    props.slider.ForegroundColor=[0 0 0];
    props.slider.Units='Char';
    props.listbox=font;
    props.listbox.FontName='interface user';
    props.listbox.FontSize=10;
    props.listbox.Style='listbox';    
    props.listbox.BackgroundColor=[0.9297 0.8711 0.7969];
    props.listbox.ForegroundColor=[0 0 0];
    props.listbox.Units='Char';
    props.listbox.HorizontalAlignment='left';
    props.button=font;
    props.button.Style='pushbutton';    
    props.button.BackgroundColor=props.window.Color;
    props.button.ForegroundColor=[0 0 0];
    props.button.Units='Char';
    props.button.BusyAction='cancel';
    props.buttonActive=font;
    props.buttonActive.Style='pushbutton';    
    props.buttonActive.BackgroundColor=.8*props.window.Color;
    props.buttonActive.ForegroundColor=[0.5 0 0];
    props.buttonActive.Units='Char';
    props.buttonActive.BusyAction='cancel';
    props.radio=font;
    props.radio.Style='radiobutton';    
    props.radio.BackgroundColor=props.window.Color;
    props.radio.ForegroundColor=[0 0 0];
    props.radio.Units='Char';
    props.edit=font;
    props.edit.Style='edit';    
    props.edit.BackgroundColor=props.listbox.BackgroundColor;
    props.edit.ForegroundColor=[0 0 0];
    props.edit.HorizontalAlignment='right';
    props.edit.Units='Char';
    props.checkbox=font;
    props.checkbox.Style='checkbox';
    props.checkbox.HorizontalAlignment='left';
    props.checkbox.BackgroundColor=props.window.Color;
    props.checkbox.ForegroundColor=[0 0 0];
    props.checkbox.Units='Char';
    props.popup=font;
    props.popup.Style='popup';
    props.popup.HorizontalAlignment='center';
    props.popup.BackgroundColor=props.window.Color;
    props.popup.ForegroundColor=[0 0 0];
    props.popup.Units='Char';
    props.text=font;
    props.text.Style='text';
    props.text.HorizontalAlignment='left';
    props.text.BackgroundColor=props.window.Color;
    props.text.ForegroundColor=[0 0 0];
    props.text.Units='Char';
    props.normaltext=font;
    props.normaltext.HorizontalAlignment='left';
    props.normaltext.Color=[0 0 0];
    props.titletext=font;
    props.titletext.FontWeight='bold';
    props.titletext.HorizontalAlignment='center';
    props.titletext.Color=[0 0 0];
    props.msgbox=font;
    props.msgbox.FontSize=get(0,'factoryUIControlFontSize');
    props.msgbox.BackgroundColor=props.window.Color;
    props.msgbox.ForegroundColor=[0 0 0];
    props.msgboxwin.Color=props.msgbox.BackgroundColor;
    v=version;
    if str2num(v(findstr(v,'(R')+2:findstr(v,')')-1))<=12;
      props.msgboxwin.Colormap=[[0 0 0];[props.msgbox.BackgroundColor]];
    end
    props.logo.Visible='off';
    props.logo.Units='Char';
    props.logo.Color=props.window.Color;
    props.axes.Color=[.97 .97 1];
    props.axes.Units='Char';
    props.axes.Box='on';
    props.axes.yLimMode='auto';
    props.axes.Layer='top';
    props.patch.FaceColor=[0.85 0.85 1];
    props.patch.EdgeColor=props.patch.FaceColor;
    props.glass.Color=[1 0 0];
    props.glass.HorizontalAlignment='center';
    props.glass.FontName='Courier';
    props.glass.FontWeight='Bold';
    props.glass.Visible='off';
%    props.glasspointer='addzero';
%    props.glasspointer='help';
    props.glasspointer='glass';
    props.textindex=font;
    props.textindex.Color=[0 0 0];
    props.textindex.Position=[0 1 0];
    props.textindex.Visible='Off';
    props.textindex.FontWeight='Bold';
    props.textindex.FontSize=12;
    props.textindex.VerticalAlignment='bottom';
    props.textindex.Erase='xor';
    props.axesindex.Color=[1 1 .6];
    props.axesindex.Units='Pixel';
    props.axesindex.XTickLabel=[];
    props.axesindex.YTickLabel=[];
    props.axesindex.XTick=[];
    props.axesindex.YTick=[];
    props.axesindex.Box='on';
    props.axesindex.Position=[0 0 0.0001 0.0001];
    props.index.xoffset=4;
    props.index.yoffset=12;
    props.marker.Marker='.';
    props.marker.MarkerSize=6;
    props.marker.LineStyle='none';
    props.marker.Color=[0 0 .8];
    props.line.Marker='none';
    props.line.LineWidth=.7;
    props.line.LineStyle='-';
    props.line.Color=[0 0 .8];
set(0,'defaultUIControlBackgroundColor',props.msgbox.BackgroundColor,...
      'defaultUIControlForegroundColor',props.msgbox.ForegroundColor,...
      'defaultTextColor',props.text.ForegroundColor)
