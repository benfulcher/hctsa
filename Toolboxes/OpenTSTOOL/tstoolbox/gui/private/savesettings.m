function savesettings(handles)

global TSTOOLdatapath;


lbox=get(handles.lboxhandle,'UserData');
loadh=get(handles.loadhandle,'UserData');
lall=get(handles.lallhandle,'UserData');
recopt=get(handles.recopthandle,'UserData');
plotmode1=get(handles.plotmodehandle(1),'Checked');
plotmode2=get(handles.plotmodehandle(2),'Checked');
defwindow=get(handles.defwindowhandle,'UserData');
surrogate=get(handles.surrogate,'UserData');
trevtc3=get(handles.trevtc3,'UserData');
position=get(handles.fighandle,'Position');


save(fullfile(TSTOOLdatapath,'tstool.mat'),'lbox','loadh','lall','recopt','plotmode1','plotmode2','defwindow',...
    'surrogate','trevtc3','position');





