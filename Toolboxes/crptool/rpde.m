function out=rpde(varargin)
% RPDE   Computes the recurrence time entropy.
%    Y=RPDE(X [,Y] [,param1,param2,...]) 
%    Calculates the normalised entropy Y of the
%    recurrence time distribution of time series X,
%    also known as recurrence period density entropy (RPDE).
%
%    Note: In contrast to the calculation of RPDE here,
%          in CRQA a Theiler window is applied to the RP
%          by default, resulting in different RPDE values. 
%          For comparison, you should ensure that the 
%          Theiler window in CRQA is set to 0.
%
%    Examples: a = sin(0:.1:80);
%              b = sin(0:.1:80) + .1 * randn(1,801);
%              rpde(a,3,15,.1)
%              rpde(b,3,15,.1)
%
%    See also CRQA, TT.
%
%    References: 
%    Little, M., McSharry, P., Roberts, S., Costello, D., Moroz, I.:
%    Exploiting Nonlinear Recurrence and Fractal Scaling Properties 
%    for Voice Disorder Detection, Biomed. Eng. Online, 6, 2007.

% Copyright (c) 2010-
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% $Date: 2011/02/24 11:36:04 $
% $Revision: 5.3 $
%
% $Log: rpde.m,v $
% Revision 5.3  2011/02/24 11:36:04  marwan
% bugfix in lmin check
% adding note on Theiler window (regarding the calculation in crqa.m)
%
% Revision 5.2  2010/06/30 12:02:31  marwan
% Help text modified
%
% Revision 5.1  2010/06/21 10:56:20  marwan
% initial submission
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% programme properties

global props

init_properties

lmin=1; % minimum length of white vertical lines to be considered
w=[]; method='max'; method_n=1; t=1; m=1; e=.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check the input

error(nargchk(1,10,nargin));
if nargout>1, error('Too many output arguments'), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check and read the input

   varargin{11}=[];
   % transform any int to double
   intclasses = {'uint8';'uint16';'uint32';'uint64';'int8';'int16';'int32';'int64'};
   flagClass = [];
   for i = 1:length(intclasses)
       i_int=find(cellfun('isclass',varargin,intclasses{i}));
       if ~isempty(i_int)
           for j = 1:length(i_int)
               varargin{i_int(j)} = double(varargin{i_int(j)});
           end
           flagClass = [flagClass; i_int(:)];
       end
   end
   if ~isempty(flagClass)
       disp(['Warning: Input arguments at position [',num2str(flagClass'),'] contain integer values']);
       disp(['(now converted to double).'])
   end
   i_double=find(cellfun('isclass',varargin,'double'));
   i_char=find(cellfun('isclass',varargin,'char'));
   nogui=0;

   if nargin & isnumeric(varargin{1})

     % check the text input parameters for method, gui 
      check_meth={'ma','eu','mi','nr','rr','fa','in','om','op','di'}; 	% maxnorm, euclidean, nrmnorm,  fan, distance
      check_gui={'gui','nog','sil'};		% gui, nogui, silent
      temp_meth=0;
      temp_gui=0;
      if ~isempty(i_char)
         for i=1:length(i_char), 
            varargin{i_char(i)}(4)='0';
            temp_gui=temp_gui+strcmpi(varargin{i_char(i)}(1:3),check_gui'); 
            temp_meth=temp_meth+strcmpi(varargin{i_char(i)}(1:2),check_meth'); 
         end
         method_n=min(find(temp_meth));
         nogui=min(find(temp_gui))-1;
         for i=1:length(i_char); temp2(i,:)=varargin{i_char(i)}(1:3); end
         i_char(strmatch(check_gui(find(temp_gui)),temp2))=[];
         if isempty(nogui), nogui=0; end
         if isempty(method_n), method_n=1; end
         if nogui>2, nogui=1; end
         if method_n>length(check_meth), method0=length(check_meth); end
         method=check_meth{method_n};
      else
         nogui=0;
         if nargout
            nogui=1;
            action='compute'; 
         end
      end
      if nogui==0
        action='init';
      else
        action='compute'; 
      end
      
      % get the parameters for creating RP
      if max(size(varargin{1}))<=3
         error('To less values in data X.')
      end
      x=double(varargin{1});
      if isempty(varargin{2}) | ~isnumeric(varargin{2}), y=x; else
      y=double(varargin{2}); end
      if sum(double(diff(x(:,1))<=0)), embed_flag=0; end
   
      if (isnumeric(varargin{2}) & max(size(varargin{2}))==1) | ~isnumeric(varargin{2})
        y=x;
        if ~isempty(varargin{i_double(2)}), m=varargin{i_double(2)}(1); else m=1; end
        if ~isempty(varargin{i_double(3)}), t=varargin{i_double(3)}(1); else t=1; end
        if ~isempty(varargin{i_double(4)}), e=varargin{i_double(4)}(1); else e=.1; end
        if ~isempty(varargin{i_double(5)}), w=varargin{i_double(5)}(1); else w=varargin{i_double(5)}; end
%        if ~isempty(varargin{i_double(6)}), wstep=varargin{i_double(6)}(1); else wstep=1; end
      else
        if ~isempty(varargin{i_double(3)}), m=varargin{i_double(3)}(1); else m=1; end
        if ~isempty(varargin{i_double(4)}), t=varargin{i_double(4)}(1); else t=1; end
        if ~isempty(varargin{i_double(5)}), e=varargin{i_double(5)}(1); else e=.1; end
        if ~isempty(varargin{i_double(6)}), w=varargin{i_double(6)}(1); else w=varargin{i_double(6)}; end
%        if ~isempty(varargin{i_double(7)}), wstep=varargin{i_double(7)}(1); else wstep=1; end
      end
    else
      error('No valid arguments.')
    end

   
    Nx=length(x); Ny=length(y);
    if size(x,1)<size(x,2), x=x'; end
    if size(y,1)<size(y,2), y=y'; end
   
   if size(x,2)>=2
      xscale=x(:,1); 
      if ~isempty(find(diff(xscale)<0)), embed_flag=0;end
   else
      xscale=(1:length(x))'; 
   end
   if size(y,2)>=2
      yscale=y(:,1); 
      if ~isempty(find(diff(yscale)<0)), embed_flag=0;end
   else
       yscale=(1:length(y))';
   end
      
      if max(size(x))~=max(size(y)),
        if ~nogui, errordlg('Data must have the same length.','Check Data'), else error('Data must have the same length.'), end
      end
      if e<0, 
        e=1; 
	if ~nogui
	   warndlg('The threshold size E can not be negative and is now set to 1.','Check Data')
	   h=findobj('Tag','crqa_eps');
	   set(h(1),'String',str2num(e)) 
        else 
	   disp('The threshold size E can not be negative and is now set to 1.'), 
        end
      end
      if t<1, 
        t=1; 
	if ~nogui
	   warndlg('The delay T can not be smaller than one and is now set to 1.','Check Data')
	   h=findobj('Tag','crqa_maxLag');
	   set(h(1),'String',str2num(t)) 
	else
	   disp('The delay T can not be smaller than one and is now set to 1.')
	end
      end
      if isempty(w), w=.5*Nx; wstep=1; end
%      if w<2, 
%        w=2; 
%	if ~nogui, warndlg('The window size W exceeds the valid range.','Check Data')
%	else, disp('The window size W exceeds the valid range.'), end
%      end
      if w>Nx, 
        w=Nx; wstep=1;; 
	if ~nogui, warndlg('The window size W exceeds the valid range.','Check Data')
	else, disp('The window size W exceeds the valid range.'), end
    end
    t=round(t); m=round(m); w=round(w);% wstep=round(wstep); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute

% calculate recurrence plot (RP)
X = crp(x,m,t,e,method,'sil');

% get recurrence times (in terms of vertical white lines in the RP)
[dummy1 dummy2 w] = tt(X);
w(w<lmin) = []; % if used with default lmin (=0), 
                 % every white line will be considered, also
                 % such of length 1 (i.e. dots)

% get histogram (normalisation will be done in ENTROPY)
h = histc(w,[0:max(w)]+.5);

% calculate entropy
out = entropy(h) / log(max(w));

