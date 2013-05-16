function Answer=inputdlg(Prompt, Title, NumLines, DefAns)
%INPUTDLG Input dialog box.
%  Answer = inputdlg(Prompt) creates a modal dialog box that returns
%  user input for multiple prompts in the cell array Answer.  Prompt
%  is a cell array containing the Prompt strings.
%
%  Answer = inputdlg(Prompt,Title) specifies the Title for the dialog.
%
%  Answer = inputdlg(Prompt,Title,LineNo) specifies the number of lines
%  for each answer in LineNo.  LineNo may be a constant value or a 
%  vector having one element per Prompt.
%
%  Answer = inputdlg(Prompt,Title,LineNo,DefAns) specifies the default
%  answer to display for each Prompt.  DefAns must contain the same
%  number of elements as Prompt and must be a cell array.
%
%  Example:
%  prompt={'Enter the matrix size:','Enter the colormap name:'};
%  def={'20','hsv'};
%  title='Input for Peaks function';
%  lineNo=1;
%  answer=inputdlg(prompt,title,lineNo,def);
%
%  See also TEXTWRAP.

%  Loren Dean   May 24, 1995.
%  Copyright (c) 1984-97 by The MathWorks, Inc.
%  $Revision: 1.1.1.1 $

%%%%%%%%%%%%%%%%%%%%%
%%% General Info. %%%
%%%%%%%%%%%%%%%%%%%%%
Black      =[0       0        0      ]/255;
LightGray  =[192     192      192    ]/255;
LightGray2 =[160     160      164    ]/255;
MediumGray =[128     128      128    ]/255;
White      =[255     255      255    ]/255;

%%%%%%%%%%%%%%%%%%%%
%%% Nargin Check %%%
%%%%%%%%%%%%%%%%%%%%
if nargout~=1,error('Wrong number of output arguents for INPUTDLG');end
if nargin<1,error('Too few arguments for INPUTDLG');end

if nargin==1,Title=' ';end
if nargin<=2, NumLines=1;end

if ~iscell(Prompt),
  Prompt={Prompt};
end

NumQuest=prod(size(Prompt));    

if nargin<=3, 
  DefAns=cell(NumQuest,1);
  for lp=1:NumQuest, DefAns{lp}=''; end
end
if nargin>4,error('Too many input arguments');end

% Backwards Compatibility
if isstr(NumLines),
  warning(['Please see the INPUTDLG help for correct input syntax.' 10 ...
           '         OKCallback no longer supported.' ]);
  NumLines=1;
end

if length(NumLines)~=NumQuest,NumLines=ones(NumQuest,1)*NumLines(1);end

if ~iscell(DefAns),
  error('Default Answer must be a cell array in INPUTDLG.');  
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% Create InputFig %%%
%%%%%%%%%%%%%%%%%%%%%%%
FigWidth=500;FigHeight=100;
FigPos(3:4)=[FigWidth FigHeight];
FigColor=get(0,'Defaultuicontrolbackgroundcolor');
InputFig=dialog(                            ...
               'Visible'         ,'off'   , ...
               'Name'            ,Title   , ...
               'Pointer'         ,'arrow' , ...
               'Units'           ,'points', ...
               'UserData'        ,''      , ...
               'Tag'             ,Title   , ...
               'HandleVisibility','on'    , ...               
               'Color'           ,FigColor, ...               
               'NextPlot'        ,'add' );
  

%%%%%%%%%%%%%%%%%%%%%
%%% Set Positions %%%
%%%%%%%%%%%%%%%%%%%%%
DefOffset=5;
SmallOffset=2;

DefBtnWidth=50;
BtnHeight=20;
BtnYOffset=DefOffset;
BtnFontSize=get(0,'FactoryUIControlFontSize');
BtnWidth=DefBtnWidth;
TxtBackClr=FigColor;
TxtForeClr=Black;

StInfo.Style              ='text'     ;
StInfo.Units              ='points'   ;   
StInfo.FontSize           =BtnFontSize;
StInfo.HorizontalAlignment='left'     ;
StInfo.BackgroundColor    =TxtBackClr ;
StInfo.ForegroundColor    =TxtForeClr ;
StInfo.HandleVisibility   ='callback' ;


EdInfo=StInfo;
EdInfo.Style='edit';
EdInfo.BackgroundColor=White;

BtnInfo=StInfo;
BtnInfo.Style='pushbutton';
BtnInfo.HorizontalAlignment='center';

% Determine # of lines for all Prompts
ExtControl=uicontrol(StInfo, ...
                     'String'   ,''         , ...    
                     'Position' ,[DefOffset                  DefOffset  ...
                                 0.95*(FigWidth-2*DefOffset) BtnHeight  ...
                                ]            , ...
                     'Visible'  ,'off'         ...
                     );
                     
WrapQuest=cell(NumQuest,1);
QuestPos=zeros(NumQuest,4);

for ExtLp=1:NumQuest,
  [WrapQuest{ExtLp},QuestPos(ExtLp,1:4)]= ...
            textwrap(ExtControl,Prompt(ExtLp),100);
end % for ExtLp
delete(ExtControl);

QuestHeight=QuestPos(:,4);

TxtHeight=QuestHeight(1)/size(WrapQuest{1,1},1);
EditHeight=TxtHeight*NumLines;
EditHeight(NumLines==1)=EditHeight(NumLines==1)+4;

FigHeight=(NumQuest+2)*DefOffset    + ...
          BtnHeight+sum(EditHeight) + ...
          sum(QuestHeight);%          + ...
                           %NumQuest*SmallOffset;

FigWidth=max(300,max(QuestPos(:,3)))+2*DefOffset;				% hier ansetzen wenn die Breite des Dialoges geaendert werden soll
FigPos=get(InputFig,'Position');

Temp=get(0,'Units');
set(0,'Units','points');
ScreenSize=get(0,'ScreenSize');
set(0,'Units',Temp);

FigPos(1)=(ScreenSize(3)-FigWidth)/2;
FigPos(2)=(ScreenSize(4)-FigHeight)/2;
FigPos(3)=FigWidth;
FigPos(4)=FigHeight;
set(InputFig,'Position',FigPos);

TxtXOffset=DefOffset;
TxtWidth=FigWidth-2*DefOffset;
TxtForeClr=Black;
TxtBackClr=get(InputFig,'Color');

QuestYOffset=zeros(NumQuest,1);
EditYOffset=zeros(NumQuest,1);
QuestYOffset(1)=FigHeight-DefOffset-QuestHeight(1);
EditYOffset(1)=QuestYOffset(1)-EditHeight(1);% -SmallOffset;

for YOffLp=2:NumQuest,
  QuestYOffset(YOffLp)=EditYOffset(YOffLp-1)-QuestHeight(YOffLp)-DefOffset;
  EditYOffset(YOffLp)=QuestYOffset(YOffLp)-EditHeight(YOffLp); %-SmallOffset;
end % for YOffLp

QuestHandle=[];
EditHandle=[];
for lp=1:NumQuest,
  QuestTag=['Prompt' num2str(lp)];
  EditTag=['Edit' num2str(lp)];
  QuestHandle(lp)=uicontrol(InputFig  ,                         ...
                           StInfo     , ...
                           'Max'      ,size(Prompt{lp},1), ...
                           'Position' ,[ TxtXOffset QuestYOffset(lp) ...
                                         TxtWidth   QuestHeight(lp)  ...
                                       ]                      , ...
                           'String'   ,WrapQuest{lp}       , ...
                           'Tag'      ,QuestTag                 ...
                           );

  EditHandle(lp)=uicontrol(InputFig   ,                       ...
                           EdInfo     , ...
                          'Max'       ,NumLines(lp)         , ...
                          'Position'  ,[ TxtXOffset EditYOffset(lp) ...
                                         TxtWidth   EditHeight(lp)  ...
                                       ]                    , ...
                          'String'    ,DefAns{lp}           , ...
                          'Tag'       ,QuestTag               ...
                          );
                                   
end % for lp

CBString='set(gcf,''UserData'',''Cancel'');uiresume';

CancelHandle=uicontrol(InputFig   ,              ...
                      BtnInfo     , ...
                      'Position'  ,[ DefOffset DefOffset ...
                                    BtnWidth  BtnHeight  ...
                                   ]           , ...
                      'String'    ,'Cancel'    , ...
                      'Callback'  ,CBString    , ...
                      'Tag'       ,'Cancel'      ...
                      );
                                   
                                   
CBString='set(gcf,''UserData'',''OK'');uiresume';

OKHandle=uicontrol(InputFig    ,              ...
                   BtnInfo     , ...
                   'Position'  ,[ FigWidth-BtnWidth-DefOffset DefOffset ...
                                  BtnWidth                    BtnHeight ...
                                ]           , ...
                  'String'     ,'OK'        , ...
                  'Callback'   ,CBString    , ...
                  'Tag'        ,'OK'          ...
                  );
    
set(InputFig ,'Visible','on');
set(findobj(InputFig),'Units','normalized','HandleVisibility','callback');
set(InputFig,'Units','points')

uiwait(InputFig);

TempHide=get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');

if any(get(0,'Children')==InputFig),
  Answer={};
  if strcmp(get(InputFig,'UserData'),'OK'),
    Answer=cell(NumQuest,1);
    for lp=1:NumQuest,
      Answer(lp)=get(EditHandle(lp),{'String'});
    end % for
  end % if strcmp
  delete(InputFig);
else,
  Answer={};
end % if any

set(0,'ShowHiddenHandles',TempHide);
