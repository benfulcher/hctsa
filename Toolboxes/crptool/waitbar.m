function fout = waitbar(x,whichbar, varargin)
%WAITBAR Display wait bar.
%   H = WAITBAR(X,'title', property, value, property, value, ...) 
%   creates and displays a waitbar of fractional length X.  The 
%   handle to the waitbar figure is returned in H.
%   X should be between 0 and 1.  Optional arguments property and 
%   value allow to set corresponding waitbar figure properties.
%
%   Property can also be an action keyword 'CreateCancelBtn', in 
%   which case a cancel button will be added to the figure, and 
%   the passed value string will be executed upon clicking on the 
%   cancel button or the close figure button.
%
%   A further helpful property is 'ShowTime', which enables to 
%   control when a field with the remaining time will appear:
%   If the passed value is 0, the remaining time will be shown,
%   if it is -1, the remaining time will never be shown, values
%   larger than 0 means, that the remaining time is shown if it
%   exceeds this value (given in seconds).
%
%   WAITBAR(X) will set the length of the bar in the most recently
%   created waitbar window to the fractional length X.
%
%   WAITBAR(X,H) will set the length of the bar in waitbar H
%   to the fractional length X.
%
%   WAITBAR(X,H,'updated title') will update the title text in
%   the waitbar figure, in addition to setting the fractional
%   length to X.
%
%   WAITBAR is typically used inside a FOR loop that performs a 
%   lengthy computation.  A sample usage is shown below:
%
%       h = waitbar(0,'Please wait...');
%       for i = 1:100,
%           % computation here %
%           waitbar(i/100,h)
%       end
%       close(h)

%   Clay M. Thompson 11-9-92
%   Vlad Kolesnikov  06-7-99
%   Copyright 1984-2000 The MathWorks, Inc.
%   $Revision: 2.6 $  $Date: 2006/02/06 13:46:17 $
%   Modified: Norbert Marwan, 2002-12-11, 2004-07-28, 2004-10-25, 2005-09-13


if nargin>=2
    if ischar(whichbar)
        type=2; %we are initializing
        name=whichbar;
    elseif isnumeric(whichbar)
        type=1; %we are updating, given a handle
        f=whichbar;
    else
        error(['Input arguments of type ' class(whichbar) ' not valid.'])
    end
elseif nargin==1
    f = findobj(allchild(0),'flat','Tag','TMWWaitbar');
    
    if isempty(f)
        type=2;
        name='Waitbar';
    else
        type=1;
        f=f(1);
    end   
else
    error('Input arguments not valid.');
end

x = max(0,min(100*x,100));

switch type
 case 1,  % waitbar(x)    update
  p = findobj(f,'Type','patch');
  pr = findobj(f,'Type','text'); pr=pr(1);
  if isempty(f) | isempty(p) | isempty(pr)
      error('Couldn''t find waitbar handles.'); 
  end

  lastUpdate = get(f,'UserData');
  dx = x;
  dt = etime(clock,lastUpdate.time);
  if dx > 0
      secRemain = dt/dx * (100 - x);
  else
      secRemain = 0;
  end
  
  timeRemain = [];
  if (secRemain > lastUpdate.showTime | lastUpdate.flag) & lastUpdate.showTime >= 0
      tx = datestr(datenum(0, 0, 0, 0, 0, secRemain),13);
      if secRemain > 86400
        tx = [datestr(datenum(0, 0, 0, 0, 0, secRemain),7),'d:',datestr(datenum(0, 0, 0, 0, 0, secRemain),13)];
      end
      timeRemain = [' (',tx,')'];
      lastUpdate.flag = 1;
  end
%  lastUpdate.time = clock;
  lastUpdate.x = x;
  set(f,'UserData',lastUpdate)

  set(pr,'String',[sprintf('%3.0f',round(x)),'%',timeRemain])
  xpatch = get(p,'XData');
  xpatch = [0 x x 0];
  set(p,'XData',xpatch)

      
            
  if nargin>2,
      % Update waitbar title:
      hAxes = findobj(f,'type','axes');
      hTitle = get(hAxes,'title');
      set(hTitle,'string',varargin{1});
  end
  
 case 2,  % waitbar(x,name)  initialize
  vertMargin = 0;
  if nargin > 2,
      % we have optional arguments: property-value pairs
      if rem (nargin, 2 ) ~= 0
          error( 'Optional initialization arguments must be passed in pairs' );
      end
  end
  
  oldRootUnits = get(0,'Units');

  set(0, 'Units', 'points');
  screenSize = get(0,'ScreenSize');
  
  axFontSize=get(0,'FactoryAxesFontSize');
  
  pointsPerPixel = 72/get(0,'ScreenPixelsPerInch');
  
  width = 360 * pointsPerPixel;
  height = 75 * pointsPerPixel;
  pos = [screenSize(3)/2-width/2 screenSize(4)/2-height/2 width height];
  
  lastUpdate.time = clock;
  lastUpdate.x = x;
  lastUpdate.flag = 0;
  showTime = 60;
  for i = 1:length(varargin)-1
    if strcmpi(varargin{i},'showTime')
        showTime = varargin{i+1};
        varargin(i:i+1) = [];
        break;
    end
  end
  
  lastUpdate.showTime = showTime;

  f = figure(...
      'Units', 'points', ...
      'DoubleBuffer', 'on', ... 
      'BusyAction', 'queue', ...
      'Position', pos, ...
      'Resize','off', ...
      'CreateFcn','', ...
      'NumberTitle','off', ...
      'IntegerHandle','off', ...
      'MenuBar', 'none', ...
      'Tag','TMWWaitbar',...
      'Interruptible', 'off', ...
      'UserData',lastUpdate, ...
      'Visible','off');
  
  %%%%%%%%%%%%%%%%%%%%%
  % set figure properties as passed to the fcn
  % pay special attention to the 'cancel' request
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if nargin > 2,
      propList = varargin(1:2:end);
      valueList = varargin(2:2:end);
      cancelBtnCreated = 0;
      for ii = 1:length( propList )
          try
              if strcmpi(propList{ii}, 'createcancelbtn' ) & ~cancelBtnCreated
                  cancelBtnHeight = 23 * pointsPerPixel;
                  cancelBtnWidth = 60 * pointsPerPixel;
                  newPos = pos;
                  vertMargin = vertMargin + cancelBtnHeight;
                  newPos(4) = newPos(4)+vertMargin;
                  callbackFcn = [valueList{ii}];
                  set( f, 'Position', newPos, 'CloseRequestFcn', callbackFcn );
                  cancelButt = uicontrol('Parent',f, ...
                                         'Units','points', ...
                                         'Callback',callbackFcn, ...
                                         'ButtonDownFcn', callbackFcn, ...
                                         'Enable','on', ...
                                         'Interruptible','off', ...
                                         'Position', [pos(3)-cancelBtnWidth*1.4, 7,  ...
                    cancelBtnWidth, cancelBtnHeight], ...
                                         'String','Cancel', ...
                                         'Tag','TMWWaitbarCancelButton');
                  cancelBtnCreated = 1;
              else
                  % simply set the prop/value pair of the figure
                  set( f, propList{ii}, valueList{ii});
              end
          catch
              disp ( ['Warning: could not set property ''' propList{ii} ''' with value ''' num2str(valueList{ii}) '''' ] );
          end
      end
  end  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
  colormap([]);
  
  axNorm=[.05 .3 .9 .2];
  axNorm=[.05 .3 .9 .28];
  axPos=axNorm.*[pos(3:4),pos(3:4)] + [0 vertMargin 0 0];
  
  h = axes('XLim',[0 100],...
           'YLim',[0 1],...
           'Box','on', ...
           'Units','Points',...
           'FontSize', axFontSize,...
           'Position',axPos,...
           'XTickMode','manual',...
           'YTickMode','manual',...
           'XTick',[],...
           'YTick',[],...
           'XTickLabelMode','manual',...
           'XTickLabel',[],...
           'YTickLabelMode','manual',...
           'YTickLabel',[]);

  tHandle=title(name);
  tHandle=get(h,'title');
  oldTitleUnits=get(tHandle,'Units');
  set(tHandle,...
      'Units',      'points',...
      'String',     name);
  
  tExtent=get(tHandle,'Extent');
  set(tHandle,'Units',oldTitleUnits);
  
  titleHeight=tExtent(4)+axPos(2)+axPos(4)+5;
  if titleHeight>pos(4)
      pos(4)=titleHeight;
      pos(2)=screenSize(4)/2-pos(4)/2;
      figPosDirty=true;
  else
      figPosDirty=false;
  end
  
  if tExtent(3)>pos(3)*1.10;
      pos(3)=min(tExtent(3)*1.10,screenSize(3));
      pos(1)=screenSize(3)/2-pos(3)/2;
      
      axPos([1,3])=axNorm([1,3])*pos(3);
      set(h,'Position',axPos);
      
      figPosDirty=true;
  end
  
  if figPosDirty
      set(f,'Position',pos);
  end

  xpatch = [0 x x 0];
  ypatch = [0 0 1 1];
   color = [0.7 0.1 0.1];

  p = patch(xpatch,ypatch,color,'EdgeColor','none','EraseMode','normal');
  

  prHandle=text(50,.35,'0%','Color',color,'FontWeight','bold','EraseMode','xor','HorizontalAlignment','center','VerticalAlignment','middle');

  
  set(f,'HandleVisibility','callback','visible','on');
  set(h,'Layer','top')
  set(prHandle,'String',[sprintf('%3.0f',round(x)),'%'])
  
  set(0, 'Units', oldRootUnits);
end  % case
drawnow;

if nargout==1,
    fout = f;
end
