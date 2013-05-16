function loadhandles(handles)

global TSTOOLdatapath;

load(fullfile(TSTOOLdatapath,'tstool.mat'));

set(handles.lboxhandle,'UserData',lbox);
set(handles.loadhandle,'UserData',loadh);
set(handles.lallhandle,'UserData',lall);
set(handles.recopthandle,'UserData',recopt);
set(handles.plotmodehandle(1),'Checked',plotmode1);
set(handles.plotmodehandle(2),'Checked',plotmode2);
set(handles.defwindowhandle,'UserData',defwindow);
set(handles.surrogate,'UserData',surrogate);
set(handles.trevtc3,'UserData',trevtc3);
set(handles.fighandle,'Position',position);
