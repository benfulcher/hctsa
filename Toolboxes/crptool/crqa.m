function xout = crqa(varargin)
%CRQA   Computes and plots the CRQA measures.
%    Y=CRQA(X [,Y] [,param1,param2,...]) 
%    Recurrence quantification analysis of the cross recurrence
%    plot of the vectors X and Y. 
%
%    The input vectors can be multi-column vectors, where
%    each column will be used as a component of the 
%    phase-space vector. However, if the first column is
%    monotonically increasing, it will be used as an
%    time scale for plotting.
%
%    Y=CRQA(X,M,T,E,W,WS,LMIN,VMIN,TW) computes the 
%    recurrence quantification analysis of the recurrence
%    plot of X by using the dimension M, delay T, the
%    size of neighbourhood E, the window size W and 
%    a window shifting value of WS. LMIN and VMIN 
%    specify the minimal length of diagonal and vertical 
%    line structures (default is 2) and TW specifies the
%    Theiler window (default is 1).
%
%    CRQA(...) without any output arguments opens a
%    GUI for interactively control the CRQA. If an
%    output is specified with using the option 'gui',
%    then the output will contain the figure handle.
% 
%    Parameters: 
%    Dimension M, delay T, the size of
%    neighbourhood E, the window size W and the shift
%    value WS are the first five numbers after the data 
%    series; if W is empty, the whole plot will be calculated.
%    The next two optional numeric parameters LMIN and VMIN
%    specify the minimal length of line structures.
%
%    As the last numeric parameter, the size of the Theiler 
%    window TW can be specified (default is 1). This window 
%    excludes the recurrence points parallel to the main 
%    diagonal from the analysis. The application of the
%    Theiler window is useful only for recurrence plots. In
%    cross recurrence plots, the size of the Theiler window will
%    be set automatically to zero.
%
%    Further parameters can be used to switch between various 
%    methods of finding the neighbours of the phasespace 
%    trajectory, to suppress the normalization of the data 
%    and to suppress the GUI (useful in order to use this 
%    programme by other programmes). The minimal length of 
%    diagonal and vertical structures can be setted only in 
%    the GUI.
%
%    Methods of finding the neighbours.
%      maxnorm     - Maximum norm.
%      euclidean   - Euclidean norm.
%      minnorm     - Minimum norm.
%      nrmnorm     - Euclidean norm between normalized vectors
%                    (all vectors have the length one).
%      rr          - Maximum norm, fixed recurrence rate.
%      fan         - Fixed amount of nearest neighbours.
%      inter       - Interdependent neighbours.
%      omatrix     - Order matrix.
%      opattern    - Order patterns recurrence plot.
%
%    Normalization of the data series.
%      normalize   - Normalization of the data.
%      nonormalize - No normalization of the data.
%
%    Suppressing the GUI.
%      gui         - Creates the GUI.
%      nogui       - Suppresses the GUI.
%      silent      - Suppresses all output.
%
%    Parameters not needed to be specified.
%
%    Output:
%      Y(:, 1) = RR     (recurrence rate)
%      Y(:, 2) = DET    (determinism)
%      Y(:, 3) = <L>    (mean diagonal line length)
%      Y(:, 4) = Lmax   (maximal diagonal line length)
%      Y(:, 5) = ENTR   (entropy of the diagonal line lengths)
%      Y(:, 6) = LAM    (laminarity)
%      Y(:, 7) = TT     (trapping time)
%      Y(:, 8) = Vmax   (maximal vertical line length)
%      Y(:, 9) = T1     (recurrence time of 1st type)
%      Y(:,10) = T2     (recurrence time of 2nd type)
%      Y(:,11) = RTE    (recurrence time entropy, i.e., RPDE)
%      Y(:,12) = Clust  (clustering coefficient)
%      Y(:,13) = Trans  (transitivity)
%
%    The window of length w is applied on the data and not on the RP, 
%    i.e. the RP will have smaller size than the window, thus w-(m-1)*tau. 
%    If we consider the data window at time i ... i+w, the corresponding RQA 
%    measures are assigned to time i. Therefore, if you see a beginning 
%    of a transition in the plot of the RQA measures at time i, this 
%    transition will probably happen at time i+w-(m-1)*tau. 
%
%    Warning:
%    The RQA measures may differ from those of the RQA programmes by
%    Charles Webber Jr. For compatibility use a Theiler window of
%    size one and ensure that the data are normalized before by the
%    same distance which is used in the RQA programmes; e.g. normalize
%    with the maximal phase space diameter, which can be estimated 
%    with the programme PSS:
%
%      RQA = crqa(100*x/pss(x,dim,lag,'euclidean'),...
%               dim,lag,e,[],[],l_min,v_min,1,...
%               'euclidean','nonormalize','silent')
%
%    Examples: a = randn(300,1);
%              crqa(a,1,1,.2,40,2,'euc')
%
%              N = 500; w = 40; ws = 2;
%              a = 3.4:.6/(N-1):4;
%              b = .5; for i = 2:N, b(i) = a(i)*b(i-1)*(1-b(i-1)); end
%              y = crqa(b,3,2,.1,w,ws);
%              subplot(2,1,1), plot(a,b,'.','markersize',.1)
%              title('logistic map'), axis([3.4 4 0 1])
%              subplot(2,1,2), plot(a(1:ws:N-w),y(1:ws:N-w,1))
%              ylabel('recurrence rate'), axis([3.4 4 0 1])
%
%      
%    See also CRQAD, CRQAD_BIG, CRP, CRP2, CRP_BIG, DL, TT, PSS.
%
%    References: 
%    Marwan, N., Romano, M. C., Thiel, M., Kurths, J.: 
%    Recurrence Plots for the Analysis of Complex Systems, 
%    Phys. Rep., 438, 2007.
%
%    Little, M., McSharry, P., Roberts, S., Costello, D., Moroz, I.:
%    Exploiting Nonlinear Recurrence and Fractal Scaling Properties 
%    for Voice Disorder Detection, Biomed. Eng. Online, 6, 2007.
%
%    Boccaletti, S., Latora, V., Moreno, Y., Chavez, M., Hwang, D.-U.:
%    Complex networks: Structures and dynamics,
%    Phys. Rep., 424, 2006.
%
%    Marwan, N., Donges, J. F., Zou, Y., Donner, R. V., Kurths, J.: 
%    Complex network approach for recurrence analysis of time series, 
%    Phys. Lett. A, 373(46), 2009. 

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 1998-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2012/11/14 13:07:03 $
% $Revision: 5.45 $
%
% $Log: crqa.m,v $
% Revision 5.45  2012/11/14 13:07:03  marwan
% bugfix for RTE for very short time series
%
% Revision 5.44  2010/09/09 10:40:01  marwan
% bugfix in cross version for calculation of network measures
%
% Revision 5.43  2010/08/27 11:16:31  marwan
% fix output of network measures
%
% Revision 5.42  2010/07/07 08:26:30  marwan
% updated references
% bug in mean clustering coefficient
%
% Revision 5.41  2010/06/30 12:01:53  marwan
% RPDE, Clustering and Transitivity added
%
% Revision 5.40  2010/06/29 12:48:56  marwan
% better debugging performance if error occurs
%
% Revision 5.39  2010/01/06 08:58:59  marwan
% adding pdist as an alternative to get the RP
%
% Revision 5.38  2009/05/14 16:20:32  marwan
% automatic window length determination restored back
% further argument checking (extreme embedding) added
%
% Revision 5.37  2009/03/31 12:28:26  marwan
% automatic window length determination corrected
%
% Revision 5.36  2009/03/24 08:35:00  marwan
% corrected DET calculation when using Theiler window
%
% Revision 5.35  2008/03/05 17:43:08  marwan
% missed Vmax
%
% Revision 5.34  2007/12/20 16:23:49  marwan
% bug in crqa output resolved
%
% Revision 5.33  2007/07/18 17:18:44  marwan
% integer values in the arguments supported
%
% Revision 5.32  2007/05/15 17:33:13  marwan
% new neighbourhood criterion: fixed RR
%
% Revision 5.31  2007/01/10 12:29:24  marwan
% bug fix: rare bug due to a Matlab bug in MenuBar property for some XWindows systems (and only if matlab is running without java)
%
% Revision 5.30  2006/11/27 11:45:40  marwan
% bug entropy vertical white lines resolved
%
% Revision 5.29  2006/10/24 14:16:16  marwan
% minor change: sigma in title line of RP shown only for normalised data
%
% Revision 5.28  2006/09/13 15:53:36  marwan
% Bug during printing solved
%
% Revision 5.27  2006/07/04 14:04:22  marwan
% order patterns for multi-column vectors without embedding
%
% Revision 5.26  2006/03/29 13:07:55  marwan
% problems regarding OPRPs and embedding resolved
%
% Revision 5.25  2006/02/14 11:45:49  marwan
% *** empty log message ***
%
% Revision 5.24  2006/02/08 13:25:12  marwan
% bug in plugin support (lmin and vmin not correct) resolved
%
% Revision 5.23  2006/02/06 14:55:41  marwan
% bug in order patterns plugin support solved
%
% Revision 5.22  2006/02/06 13:46:17  marwan
% plugin for order patterns recurrence plots supported
%
% Revision 5.21  2005/11/28 10:16:35  marwan
% && and || changed to & and |
% (seems to cause problems in Matlab 12.1)
%
% Revision 5.20  2005/09/13 12:08:12  marwan
% fix of a serious bug in the GUI
%
% Revision 5.19  2005/07/29 09:02:59  marwan
% normalization bug resolved
%
% Revision 5.18  2005/04/15 09:02:32  marwan
% minor bugfix in plugin section
%
% Revision 5.17  2005/04/08 09:54:11  marwan
% plugin added
%
% Revision 5.15  2005/04/06 13:29:12  marwan
% computation of recurrence times is now faster
%
% Revision 5.14  2005/04/06 12:59:24  marwan
% several small bug fixes which improves the precision
%
% Revision 5.13  2005/04/01 12:19:51  marwan
% bug in minimal length for diagonal and vertical lines fixed
%
% Revision 5.12  2005/03/16 13:16:12  marwan
% bug in output fixed (same time scales for all sub-plots)
%
% Revision 5.11  2005/03/16 12:23:08  marwan
% support long data series by using crp_big
%
% Revision 5.10  2005/03/16 11:23:52  marwan
% automatic detection:
% * of provided time-scale
% * if multi-column input can be used as phase-space vector
%
% Revision 5.9  2004/12/23 07:49:03  marwan
% bug in order patterns RP fixed (empty order patterns)
%
% Revision 5.8  2004/11/15 12:56:01  marwan
% bug fix in method checking
%
% Revision 5.7  2004/11/15 12:53:27  marwan
% covariance plot commented out
%
% Revision 5.6  2004/11/12 08:40:46  marwan
% order patterns recurrence plot added
%
% Revision 5.5  2004/11/10 07:04:50  marwan
% initial import
%
%
% This program is part of the new generation XXII series.
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% programme properties

global errcode props


init_properties

lmin=2;
vmin=2;
nonorm=1;
theiler_window=1;
hw=-1;
xscale = [];
yscale = [];
undocumented = 0;
argout_index=[1,2,3,4,5,6,7,8,9,10,11,12,13];
RTE = NaN; Clust = NaN; Trans = NaN;

flag_pdist = 0; % use Matlab function PDIST, but only working for single x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check the input

error(nargchk(1,13,nargin));
if nargout>1, error('Too many output arguments'), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splash the GPL

splash_gpl('crp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check and read the input

try 
errcode=1;

   set(0,'ShowHidden','on')
   varargin{14}=[];
   t=1;
   w=[];wstep=0; method='max'; method_n=1;
   time_scale_flag=1; time_scale_flag_x=1;
   nogui=0;
   
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
      check_meth={'ma','eu','mi','nr','rr','fa','in','om','op','di'}; 	% maxnorm, euclidean, nrmnorm,  fan, distance
      check_gui={'gui','nog','sil'};                % gui, nogui, silent
      check_norm={'non','nor'};                        % nonormalize, normalize
   
   if nargin & ischar(varargin{1})
        set(0,'ShowHidden','on');h=findobj('Tag','Msgbox_Check Input');delete(h)
        action=varargin{1};
        h=findobj('Tag','crqa_Fig');
        temp=get(h(1),'UserData');
        x=temp{1}; y=temp{2};
        h=findobj('Tag','crqa_axes_Data','Parent',h(1)); % get handle of plot-axes object
        h_crqa_axes_Data = h(1);
        h=findobj('Type','line','Parent',h(1)); % get handle of line-object
        xscale=get(h(1),'xdata');
        h=findobj('Tag','crqa_m');
        m=str2num(get(h(1),'String'));
        h=findobj('Tag','crqa_maxLag');
        t=str2num(get(h(1),'String'));
        h=findobj('Tag','crqa_eps');
        e=str2num(get(h(1),'String'));
        h=findobj('Tag','crqa_method');
        method={'Maximum Norm','Euclidean Norm','Minimum Norm','Normalized Norm','Maximum Norm, fixed RR','FAN','Interdependent','Order Matrix','Order Pattern','Distance Plot'};
        method=method{get(h(1),'Value')};
        method_n = get(h(1),'Value');
        nonorm=get(h(1),'UserData');
        h=findobj('Tag','crqa_lmin');
        lmin=str2num(get(h(1),'String'));
        h=findobj('Tag','crqa_vmin');
        vmin=str2num(get(h(1),'String'));
        h=findobj('Tag','crqa_theiler');
        theiler_window=str2num(get(h(1),'String'));
        h=findobj('Tag','crqa_w');
        w=str2num(get(h(1),'String'));
        h=findobj('Tag','crqa_ws');
        wstep=str2num(get(h(1),'String'));

   elseif nargin & isnumeric(varargin{1})

       % check the text input parameters for method, gui 
        temp_meth=0;
        temp_norm=0;
        temp_gui=0;
        if ~isempty(i_char)
             for i=1:length(i_char), 
                varargin{i_char(i)}(4)='0';
                temp_gui=temp_gui+strcmpi(varargin{i_char(i)}(1:3),check_gui'); 
                temp_norm=temp_norm+strcmpi(varargin{i_char(i)}(1:3),check_norm'); 
                temp_meth=temp_meth+strcmpi(varargin{i_char(i)}(1:2),check_meth'); 
             end
             method_n=min(find(temp_meth));
             nogui=min(find(temp_gui))-1;
             nonorm=min(find(temp_norm))-1;
             for i=1:length(i_char); temp2(i,:)=varargin{i_char(i)}(1:3); end
             i_char(strmatch(check_gui(find(temp_gui)),temp2))=[];
             if isempty(nonorm), nonorm=1; end
             if nonorm>1, nonorm=1; end
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
            disp('Error using ==> crqa')
            disp('To less values in data X.')
            return
        end
        x=double(varargin{1});
        if isempty(varargin{2}) | ~isnumeric(varargin{2}), y=x; else
        y=double(varargin{2}); end
  %      if sum(double(diff(x(:,1))<=0)), time_scale_flag=0; end

        if (isnumeric(varargin{2}) & max(size(varargin{2}))==1) | ~isnumeric(varargin{2})
            y=x;
            if ~isempty(varargin{i_double(2)}), m=varargin{i_double(2)}(1); else m=1; end
            if ~isempty(varargin{i_double(3)}), t=varargin{i_double(3)}(1); else t=1; end
            if ~isempty(varargin{i_double(4)}), e=varargin{i_double(4)}(1); else e=.1; end
            if ~isempty(varargin{i_double(5)}), w=varargin{i_double(5)}(1); else w=varargin{i_double(5)}; end
            if ~isempty(varargin{i_double(6)}), wstep=varargin{i_double(6)}(1); else wstep=1; end
            if ~isempty(varargin{i_double(7)}), lmin=varargin{i_double(7)}(1); end
            if ~isempty(varargin{i_double(8)}), vmin=varargin{i_double(8)}(1); end
            if ~isempty(varargin{i_double(9)}), theiler_window=varargin{i_double(9)}(1); end
        else
            if ~isempty(varargin{i_double(3)}), m=varargin{i_double(3)}(1); else m=1; end
            if ~isempty(varargin{i_double(4)}), t=varargin{i_double(4)}(1); else t=1; end
            if ~isempty(varargin{i_double(5)}), e=varargin{i_double(5)}(1); else e=.1; end
            if ~isempty(varargin{i_double(6)}), w=varargin{i_double(6)}(1); else w=varargin{i_double(6)}; end
            if ~isempty(varargin{i_double(7)}), wstep=varargin{i_double(7)}(1); else wstep=1; end
            if ~isempty(varargin{i_double(8)}), lmin=varargin{i_double(8)}(1); end
            if ~isempty(varargin{i_double(9)}), vmin=varargin{i_double(9)}(1); end
            if ~isempty(varargin{i_double(10)}), theiler_window=varargin{i_double(10)}(1); end
        end
   else
        disp('Error using ==> crqa')
        disp('No valid arguments.')
        return
   end
   if method==8 & m > 1, 
       m=1; 
       disp('Warning: For order matrix a dimension of one is used.')
   end
   
   
   Nx=length(x); Ny=length(y);
   if size(x,1)<size(x,2), x=x'; end
   if size(y,1)<size(y,2), y=y'; end


   if (strcmpi(method,'Order Pattern') | strcmpi(method,'op')) & m == 1 & size(x,2) < 2
       m=2; 
       h = findobj('Tag','crqa_m');
       if ~isempty(h)
           errordlg(['For order patterns recurrence plots the',10,'dimension must be larger than one.'],'Dimension too small')
           set(h(1),'String',num2str(m))
       else
           disp(['Warning: For order patterns recurrence plots the dimension must',10,...
                'be larger than one. ',...
                'Embedding dimension is set to ',num2str(m),'.'])
       end
   end

   if size(x,2)>=2
      xscale=x(:,1); 
      if ~isempty(find(diff(xscale)<0)) % multi-column data vector, each column used as vector component
          time_scale_flag=0; 
          time_scale_flag_x=0; 
          xscale=(1:length(x))'; 
      end
   else
      xscale=(1:length(x))'; 
      time_scale_flag=0; 
      time_scale_flag_x=0; 
   end

   if size(y,2)>=2
      yscale=y(:,1); 
      time_scale_flag=1; 
      if ~isempty(find(diff(yscale)<0))
          time_scale_flag=0; 
          yscale=(1:length(y))';
      end
      if time_scale_flag & ~time_scale_flag_x % if time-scale given in y, but not in x -> error
          if ~nogui
             errordlg(['A time-scale for the second data is series given, but not for the first!',10,'(The time-scale has to be inlcuded as the first colummn of the first data vector.)'],'Check Data')
             waitforbuttonpress
             return
          else
             error(['A time-scale for the second data series is given, but not for the first!',10,'(The time-scale has to be inlcuded as the first colummn of the first data vector.)'])
          end
      end
   else
       yscale=(1:length(y))';
       time_scale_flag=0; 
   end
   if time_scale_flag_x & ~time_scale_flag % if time-scale given in x, but not in y
       if length(x) ~= length(y)
          if ~nogui
             errordlg(['If you are using the time-scale given by the first vector also for',10,'the second vector, both vectors should have the same size!'],'Check Data')
             waitforbuttonpress
             return
          else
             error(['If you are using the time-scale given by the first vector also for',10,'the second vector, both vectors should have the same size!'])
          end
       end
       y = [xscale, y];
       yscale = xscale;
       time_scale_flag=1; 
   end

   if max(size(x))~=max(size(y)),
        if ~nogui, errordlg('Data must have the same length.','Check Data'), waitforbuttonpress, return, else error('Data must have the same length.'), end
   end
   if e<0, 
       e=1; 
       if ~nogui
          warndlg('The threshold size E can not be negative and is now set to 1.','Check Data')
          waitforbuttonpress
          h=findobj('Tag','crqa_eps');
          if ~isempty(h), set(h(1),'String',num2str(e)), end
       else 
          disp('The threshold size E can not be negative and is now set to 1.'), 
       end
   end
   if t<1, 
     t=1; 
     if ~nogui
        warndlg('The delay T can not be smaller than one and is now set to 1.','Check Data')
        waitforbuttonpress
        h=findobj('Tag','crqa_maxLag');
        if ~isempty(h), set(h(1),'String',num2str(t)), end
     else
        disp('The delay T can not be smaller than one and is now set to 1.')
     end
   end
   if m < 1; m = 1; end


   if (Nx - (m-1) * t) < 5
        disp('Error using ==> crqa')
        disp('Too less data.')
        xout = [];
        return
   end


%   if isempty(w), w=Nx-(m-1)*t; wstep=1; end
   if isempty(w), w=Nx; wstep=1; end
   if w < 5+(m-1)*t, 
     w=5+(m-1)*t;
     if w > Nx, w = Nx; end
     if ~nogui, warndlg('The window size W falls below the valid range.','Check Data')
        waitforbuttonpress
        h=findobj('Tag','crqa_w');
        if ~isempty(h), set(h(1),'String',num2str(w)), end
     else, disp('The window size W falls below the valid range.'), end
   end
%   if w>Nx-(m-1)*t, 
   if w>Nx, 
%     w=Nx-(m-1)*t; wstep=1;; 
     w=Nx; wstep=1;; 
     if ~nogui, warndlg('The window size W exceeds the valid range.','Check Data')
        waitforbuttonpress
        h=findobj('Tag','crqa_w');
        if ~isempty(h), set(h(1),'String',num2str(w)), end
     else, disp('The window size W exceeds the valid range.'), end
   end
   if wstep<1 | wstep>Nx/3, 
     wstep=2; 
     if ~nogui, warndlg('The window shifting value WS exceeds the valid range.','Check Data')
        waitforbuttonpress
        h=findobj('Tag','crqa_ws');
        if ~isempty(h), set(h(1),'String',num2str(wstep)), end
     else, disp('The window shifting value WS exceeds the valid range.'), end
   end
   if vmin<1 | vmin>Nx, 
     vmin=2; 
     if ~nogui, warndlg('The minimal length for vertical lines is not valid.','Check Data')
        waitforbuttonpress
        h=findobj('Tag','crqa_vmin');
        if ~isempty(h), set(h(1),'String',num2str(vmin)), end
     else, disp('The minimal length for vertical lines is not valid.'), end
   end
   if lmin<1 | lmin>Nx, 
     lmin=2; 
     if ~nogui, warndlg('The minimal length for diagonal lines is not valid.','Check Input')
        waitforbuttonpress
        h=findobj('Tag','crqa_lmin');
        if ~isempty(h), set(h(1),'String',num2str(lmin)), end
     else, disp('The minimal length for diagonal lines is not valid.'), end
   end
   if theiler_window<0 | theiler_window>Nx, 
     theiler_window=1; 
     if ~nogui, warndlg('The value for the Theiler window is not valid.','Check Input')
        waitforbuttonpress
        h=findobj('Tag','crqa_theiler');
        if ~isempty(h), set(h(1),'String',num2str(theiler_window)), end
     else, disp('The value for the Theiler window is not valid.'), end
   end
   if x~=y, theiler_window=0; end
   t=round(t); m=round(m); w=round(w); wstep=round(wstep); vmin=round(vmin); lmin=round(lmin); theiler_window=round(theiler_window);
   

   % normalize data if necessary
   if nonorm
       x = (x - repmat(mean(x),length(x),1)) ./ repmat(std(x),length(x),1);
       y = (y - repmat(mean(y),length(y),1)) ./ repmat(std(y),length(y),1);
   end
      

switch(action)

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% create gui

case 'init'

  errcode=3;
       
  oldunit=get(0,'unit');
  set(0,'Unit','char')
  scr=get(0,'ScreenSize'); 
  set(0,'Unit',oldunit)
  
  h8=figure(props.window,...                                % Plot Figure
            'Tag','crqa_Fig',...
            'MenuBar','Figure',...
            'Position',[(scr(3)-150)/2 scr(4)-50 150.0000 44],...
            'PaperType','a4',...
            'PaperPosition',[0.25 0.25 7.7677 11.193],...
            'PaperOrientation','portrait',...
            'MenuBar','none', ...
            'UserData',{x,y,time_scale_flag},...
            'Name','Cross Recurrence Quantification Analysis');
  
  set(0,'showhidden','on')
  h=findobj('Label','&Help','Type','uimenu');
  if isempty(h)
    h=uimenu('Label','&Help');
    h2=uimenu('Parent',h(1),'Label','&Help Cross Recurrence Quantification Analysis','Callback','helpwin crqa');
  else
    h1=flipud(get(h(1),'Children'));
    set(h1(1),'Separator','on')
    h2=uimenu('Parent',h(1),'Label','&Help Cross Recurrence Quantification Analysis','Callback','helpwin crqa');
    copyobj(h1,h(1))
    delete(h1)
  end
  set(0,'showhidden','off')


  h=axes(props.axes,...
            'Position',[89+30 24.8+14.5 6.8 3.5]);    
  logo=load('logo');
  h2=imagesc([logo.logo fliplr(logo.logo)]);
  set(h2,'Tag','uniLogo')
  set(h,props.logo,'Tag','axes_logo')
  h=uicontrol(props.text,...
            'Tag','text_logo',...
            'String','Uni Potsdam',...
            'Position',[97+30 24.2143+14.5 22  3.5714]);    
  h2=textwrap(h,{[char(169),' AGNLD'],'University of Potsdam','1998-2007'});
  set(h,'String',h2)


%%%%%%%%%%% plots

  axes_height = 5.8;
  axes_base = 3;
  axes_hoffset = 2.5;
  h=axes(props.axes,...
            'Tag','crqa_axes_Data',...
                'Box','On',...
            'Position',[11.2017  axes_base+4*(axes_height+axes_hoffset)    28.1785+15    axes_height]);    

  if time_scale_flag
      plot(xscale,x(:,2:end),'color',props.line.Color)
      if ~all(x(:)==y(:)) hold on; plot(yscale,y(:,2:end),'r'), end
  else
      plot(xscale,x(:,1:end),'color',props.line.Color)
      if ~all(x(:)==y(:)) hold on; plot(yscale,y(:,1),'r'), end
  end
  set(h,'Tag','crqa_axes_Data','color',props.axes.Color)
  ylabel('Data')

  h=axes(props.axes,...
            'Tag','crqa_axes_Var',...
            'Box','On',...
            'Position',[49.8023+15   axes_base+4*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  if time_scale_flag
     x_var = winplot(x(:,2:end),w,wstep,2);
  else
     x_var = winplot(x(:,1),w,wstep,2);
  end
  if size(x_var,1) > 2
    h2=stairs(xscale(round(x_var(:,1))),x_var(:,2));
    set(h2,'color',props.line.Color)
  else
    cla
    text(0.5,0.5,sprintf('%6.4f',x_var(1,2)),'FontWeight','bold','HorizontalAlign','Center')
  end
  if ~all(x(:)==y(:)) 
    hold on; 
    if time_scale_flag
        y_var = winplot(y(:,2:end),w,wstep,2); 
    else
        y_var = winplot(y(:,1),w,wstep,2); 
    end
    if size(y_var,1) > 2
      stairs(yscale(round(y_var(:,1))),y_var(:,2),'r')
    else
      cla
      text(0.5,0.5,sprintf('%6.5f, %6.5f',x_var(2,2),y_var(1,2)),'FontWeight','bold','HorizontalAlign','Center')
    end
    h3=axes('Units','Char',...
              'Pos',get(h,'Pos'),...
              'XLim',get(h,'XLim'),...
              'YLim',get(h,'YLim'),...
              'YAxisLocation','right',...
              'Color','none',...
              'Tag','crqa_axes_CoVar',...
              'visible','off');
  end
  set(gcf,'CurrentAxes',h);
  ylabel('Variance')
  set(h,'Tag','crqa_axes_Var','color',props.axes.Color)


  h=axes(props.axes,...
            'Tag','crqa_axes_RR',...
            'Box','On',...
            'Position',[11.2017   axes_base+3*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('RR')

  h=axes(props.axes,...
            'Tag','crqa_axes_DET',...
            'Box','On',...
            'Position',[49.8023+15   axes_base+3*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('DET')

  h=axes(props.axes,...
            'Tag','crqa_axes_L',...
            'Box','On',...
            'Position',[11.2017   axes_base+2*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('L')

  h=axes(props.axes,...
            'Tag','crqa_axes_ENTR',...
            'Box','On',...
            'Position',[49.8023+15   axes_base+2*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('ENTR')

  h=axes(props.axes,...
            'Tag','crqa_axes_LAM',...
            'Box','On',...
            'Position',[11.2017    axes_base+1*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('LAM')

  h=axes(props.axes,...
            'Tag','crqa_axes_TT',...
            'Box','On',...
            'Position',[49.8023+15    axes_base+1*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('TT')

  h=axes(props.axes,...
            'Tag','crqa_axes_RTE',...
            'Box','On',...
            'Position',[11.2017    axes_base+0*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('T_1')

  h=axes(props.axes,...
            'Tag','crqa_axes_T2',...
            'Box','On',...
            'Position',[49.8023+15    axes_base+0*(axes_height+axes_hoffset)   28.1785+15    axes_height]);    
  ylabel('T_2')

%%%%%%%%%%% embedding
  h=uicontrol(props.frame,...
            'Tag','frame',...
            'Position',[86+30 29.+2.8 29 5.7]);    
 
  h=uicontrol(props.text,...
            'Tag','text',...
            'Fontangle','italic',...
            'String','Embedding parameters',...
            'Position',[87+30 34.2+1.6 16.8333  1.5000]);    
 
  h=uicontrol(props.text,...
            'Tag','text',...
            'String','Dimension:',...
            'Position',[89+30 32+2.2 16.8333  1.5000]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_m',...
            'String',num2str(m),...
            'ToolTip','Select the embedding dimension.',...
            'Position',[104+30 32+.2+2.2 7  1.5000]);    

  h=uicontrol(props.text,...
            'Tag','text',...
            'String','Delay:',...
            'Position',[89+30 30+2.32 16.8333  1.5000]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_maxLag',...
            'String',num2str(t),...
            'ToolTip','Insert the embedding delay time.',...
            'Position',[104+30 30+.2+2.32 7  1.5000]);    

%%%%%%%%%%% neigbourhood
  h=uicontrol(props.frame,...
            'Tag','frame',...
            'Position',[86+30 21.9+2.5 29 6.7]);    
 
  h=uicontrol(props.text,...
            'Tag','text',...
            'Fontangle','italic',...
            'String','Neighbourhood',...
            'Position',[87+30 26.8+2.5 23  1.5000]);    
 
  h=uicontrol(props.popup,...
            'Tag','text',...
            'Tag','crqa_method',...
            'UserData',nonorm,...
            'Value',method_n,...
            'String','Maximum Norm|Euclidean Norm|Minimum Norm|Normalized Norm|Fixed RR|Fixed Amount|Interdependent|Order Matrix|Order Patterns|Distance Plot',...
            'Position',[89+30 24.6+.2+2.5 22  1.7]);    

  h=uicontrol(props.text,...
            'Tag','text',...
            'String','Threshold:',...
            'Position',[89+30 22.6+2.5 16.8333  1.5]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_eps',...
            'String',num2str(e),...
            'ToolTip','Insert the size of neighbourhood.',...
            'Position',[104+30 22.6+.2+2.5 7  1.5000]);    

%%%%%%%%%%% crqa parameters
  h=uicontrol(props.frame,...
            'Tag','frame',...
            'Position',[86+30 11+.5 29 12.2]);    
 
  h=uicontrol(props.text,...
            'Tag','text',...
            'Fontangle','italic',...
            'String','CRQA parameters',...
            'Position',[87+30 19.4+2.5 23  1.5000]);    

  h=uicontrol(props.text,...
            'Tag','text',...
            'String','min. Diagonal:',...
            'Position',[89+30 17.6+2.5 16.8333  1.5]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_lmin',...
            'String',num2str(lmin),...
            'ToolTip','Insert the minimal length of a diagonal line.',...
            'Position',[104+30 17.6+.2+2.5 7  1.5000]);    

  h=uicontrol(props.text,...
            'Tag','text',...
            'String','min. Vertical:',...
            'Position',[89+30 15.6+2.5 16.8333  1.5]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_vmin',...
            'String',num2str(vmin),...
            'ToolTip','Insert the minimal length of a vertical line.',...
            'Position',[104+30 15.6+.2+2.5 7  1.5000]);    

  h=uicontrol(props.text,...
            'Tag','text',...
            'String','Theiler wind.:',...
            'Position',[89+30 15.6+.5 16.8333  1.5]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_theiler',...
            'String',num2str(theiler_window),...
            'ToolTip','Insert the size for the Theiler window.',...
            'Position',[104+30 15.6+.2+.5 7  1.5000]);    
  if x~=y; set(h,'enable','off'); end

  h=uicontrol(props.text,...
            'Tag','text',...
            'String','Window size:',...
            'Position',[89+30 13.6+.5 16.8333  1.5]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_w',...
            'String',num2str(w),...
            'ToolTip','Insert the size of the sliding window.',...
            'Position',[104+30 13.6+.2+.5 7  1.5000]);    

  h=uicontrol(props.text,...
            'Tag','text',...
            'String','Window step:',...
            'Position',[89+30 11.6+.5 16.8333  1.5]);    

  h=uicontrol(props.edit,...
            'Tag','crqa_ws',...
            'String',num2str(wstep),...
            'ToolTip','Insert the step width for sliding the window.',...
            'Position',[104+30 11.6+.2+.5 7  1.5000]);    


%%%%%%%%%%% buttons

  h=uicontrol(props.frame,...
            'Tag','frame',...
            'Position',[86+30 1.3+.5 29 9]);    


  h=uicontrol(props.button,...
                'String','Store',...
                'Tag','crqa_button_store',...
              'Enable','Off',...
                'ToolTip','Stores the CRQA analysis into a variable in the workspace.',...
                'Callback','crqa store',...
                'Position',[100.5+30  7.4+.5 10.5  2.2143]);

  h=uicontrol(props.button,...
            'Tag','crqa_button_print',...
            'CallBack','crqa print',...
                'ToolTip','Prints the CRQA window.',...
            'String','Print',...
            'Position',[89+30  7.4+.5 10.5  2.2143]);    

  h=uicontrol(props.button,...
            'Tag','crqa_button_close',...
            'CallBack','crqa close',...
                'ToolTip','Closes the CRQA window.',...
            'String','Close',...
            'Position',[89+30  2.1+.5 22  2.2143]);    
  
  h=uicontrol(props.button,...
                'String','Apply',...
                'Tag','crqa_button_apply',...
                'ToolTip','Starts the computation.',...
                'Callback','crqa compute',...
                'Position',[89+30  4.75+.5 22  2.2143]);
  
  set(0,'ShowHidden','on')
  set(h8, 'HandleVis','CallBack')
  tags={'crqa_Fig';'axes_logo';'text_logo';'crqa_theiler';'frame';'text';'crqa_axes_Data';'crqa_axes_Var';'crqa_axes_CoVar';'crqa_axes_RR';'crqa_axes_DET';'crqa_axes_L';'crqa_axes_ENTR';'crqa_axes_LAM';'crqa_axes_TT';'crqa_axes_RTE';'crqa_axes_T2';'crqa_m';'crqa_maxLag';'crqa_method';'crqa_eps';'crqa_lmin';'crqa_vmin';'crqa_w';'crqa_ws';'crqa_button_store';'crqa_button_print';'crqa_button_close';'crqa_button_apply'};
  h=[];
  for i=1:length(tags); h=[h; findobj('Tag',tags{i})]; end
  set(h,'Units','Norm')
  if nargout; xout=h8; end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% close windows 

case 'close'
  errcode=101;
  set(0,props.root)
  h=findobj('Tag','crqa_Fig');
  if ~isempty(h), close(h(1)), end 
  clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% store

case 'store'
  errcode=0;
  if ~isempty(findobj('Tag','crqa_button_store'))
    h=findobj('Tag','crqa_button_store');
    h1=findobj('Tag','crqa_button_close');
    if ~isempty(h1), vname_old=get(h1(1),'UserData'); else vname_old=''; end
    if isempty(vname_old), vname_old=''; end
    vname=char(inputdlg('Choose a variable name.','Store output',1,{vname_old}));
    if isempty(vname)
      return
    else
      crqa_values=get(h(1),'UserData');
      assignin('base',vname, [crqa_values])
      warndlg(['CRQA measures have been assigned to the workspace variable ''',vname,'''.'],'Store output');
           waitforbuttonpress
      set(h1(1),'UserData',vname)
    end
  end
  set(0,props.root)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% print
  case 'print'

    errcode=91;
    h=findobj('Tag','crqa_axes_Data','Parent',gcf); h_axes.h(1)=h(1);
    h=findobj('Tag','crqa_axes_Var','Parent',gcf); h_axes.h(2)=h(1);
    h=findobj('Tag','crqa_axes_CoVar','Parent',gcf); if ~isempty(h), h_axes.h(11)=h(1); end
    h=findobj('Tag','crqa_axes_RR','Parent',gcf); h_axes.h(3)=h(1);
    h=findobj('Tag','crqa_axes_DET','Parent',gcf); h_axes.h(4)=h(1);
    h=findobj('Tag','crqa_axes_L','Parent',gcf); h_axes.h(5)=h(1);
    h=findobj('Tag','crqa_axes_ENTR','Parent',gcf); h_axes.h(6)=h(1);
    h=findobj('Tag','crqa_axes_LAM','Parent',gcf); h_axes.h(7)=h(1);
    h=findobj('Tag','crqa_axes_TT','Parent',gcf); h_axes.h(8)=h(1);
    h=findobj('Tag','crqa_axes_RTE','Parent',gcf); h_axes.h(9)=h(1);
    h=findobj('Tag','crqa_axes_T2','Parent',gcf); h_axes.h(10)=h(1);
    h=findobj('Tag','uniLogo');
    tags={'text_logo';'frame';'text';'crqa_m';'crqa_maxLag';'crqa_method';'crqa_eps';'crqa_lmin';'crqa_vmin';'crqa_theiler';'crqa_w';'crqa_ws';'crqa_button_store';'crqa_button_print';'crqa_button_close';'crqa_button_apply';};
    for i=1:length(tags); h=[h; findobj('Tag',tags{i},'Parent',gcf)]; end
    set(h,'Visible','Off')
    
    set(h_axes.h,'Units','Character');
    h_axes.old_pos=get(h_axes.h,'Position');

    axes_height = .13;
    axes_base = 0.065;
    axes_hoffset = .06;
    
    for i=2:2:10
        set(h_axes.h(i-1),  'Units','normalize','Position',[0.1300    axes_base+(5-i/2)*(axes_height+axes_hoffset)    0.3270    axes_height])
        set(h_axes.h(i),'Units','normalize','Position',[0.5780    axes_base+(5-i/2)*(axes_height+axes_hoffset)    0.3270    axes_height])
    end
    if length(h_axes.h) > 10
        set(h_axes.h(11),  'Units','normalize','Position',[0.5780    axes_base+(5-2/2)*(axes_height+axes_hoffset)    0.3270    axes_height])
    end
    h_dlg=printdlg;
    waitfor(h_dlg)

    for i=1:10, set(h_axes.h(i),'Units','Character','Position',h_axes.old_pos{i}), set(h_axes.h(i),'Units','Norm'),end
    if length(h_axes.h) > 10
        set(h_axes.h(11),'Units','Character','Position',h_axes.old_pos{11}), set(h_axes.h(11),'Units','Norm')
    end
    set(h,'Visible','On')
    set(0,props.root)







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute
case 'compute'

  errcode=11;
  if length(method)>1 & strcmpi(method(1:2),'di')
      disp('Warning: RQA from distance plot not possible!')
      return
  end
  
  if ~nogui
    h_fig=findobj('tag','crqa_Fig');
    setptr(gcf,'watch'), 
    obj=({'text';'crqa_m';'crqa_maxLag';'crqa_method';'crqa_eps';'crqa_lmin';'crqa_vmin';'crqa_theiler';'crqa_w';'crqa_ws';'crqa_button_store';'crqa_button_print';'crqa_button_close'});
    for j=1:length(obj); 
        h=findobj('Tag',obj{j},'Parent',h_fig(1)); 
        if ~isempty(h)
          set(h,'Enable','Off')
        end
    end
    h=findobj('tag','crqa_button_apply');
    set(h(1),'ToolTip','Stops the computation.','String','Stop','Callback','set(0,''ShowHidden'',''on'');h=findobj(''tag'',''crqa_button_apply'');set(h(1),''String'',''Stopped'');set(0,''ShowHidden'',''off'')')
  end

if Nx==w & wstep<2, wstep=1; Nx=w+1; end
if Nx==w, Nx=w+1; end
if nogui~=2, hw=waitbar(0,['0/',num2str(Nx-w)]);set(hw,'Name','Please Wait!');h1=get(hw,'chil');h1=get(h1,'title'); drawnow; end
if strcmpi(method,'Order Pattern') method = 'op'; end
if strcmpi(method,'Order Matrx') method = 'om'; end
if strcmpi(method,'Maximum Norm, fixed RR') method = 'rr'; end


% check if plugin exist and is executable
[plugin_exist, plugin_name, plugin_path] = is_crp_plugin;
if nogui == 1 & plugin_exist & ( method_n < 4 | method_n == 8 ) & length(x) == length(y)
    disp('(plugin used)')
end

errcode=20;

% general histograms of line structures for significance test
hist_l = []; hist_v = []; hist_w = [];

for i=1:wstep:Nx-w; 
     if ~nogui
         set(0,'ShowHidden','on')
         h=findobj('tag','crqa_button_apply','Parent',h_fig(1));
         if strcmpi(get(h(1),'string'),'stopped')
              Y(i:Nx-w,1:6)=NaN;
              break
         end
     end


     if time_scale_flag
        x_var=var(x(i:i+w-1,2:end)); y_var=var(y(i:i+w-1,2:end)); 
        temp=cov(x(i:i+w-1,2:end),y(i:i+w-1,2:end));
     else
        x_var=var(x(i:i+w-1,1)); y_var=var(y(i:i+w-1,1)); 
        temp=cov(x(i:i+w-1,1),y(i:i+w-1,1));
     end
     xy_var=temp(1,2);

          
     do_norm = {'non';'nor'};

     % if plugin exist and method is MAX, MIN or EUC
     if plugin_exist & ( method_n < 4 | method_n == 9 ) & length(x) == length(y) 

         errcode=21;
         warning off
         while 1
             tmp_xdatafile = tempname;
             if ~exist(tmp_xdatafile), break, end
         end
         while 1
             tmp_ydatafile = tempname;
             if ~exist(tmp_ydatafile), break, end
         end
         while 1
             tmp_rqadatafile = tempname;
             if ~exist(tmp_rqadatafile), break, end
         end
         
         
         x_tmp = x(i:i+w-1,:);
         y_tmp = y(i:i+w-1,:);
         
         % save data in temporary file
         save(tmp_xdatafile,'x_tmp','-ascii','-tabs');
         save(tmp_ydatafile,'y_tmp','-ascii','-tabs');

         % call extern rp programme
         m_str = {'MAX', 'EUC', 'MIN', 'NR', 'RR', 'FAN', 'IN', 'OM', 'OP', 'EUC'};

         [status ] = system([plugin_path,filesep,plugin_name,' -m ',num2str(m), ...
                                       ' -t ',num2str(t), ...
                                       ' -e ',num2str(e), ...
                                       ' -n ',m_str{method_n}, ...
                                       ' -w ',num2str(theiler_window), ...
                                       ' -l ',num2str(lmin), ...
                                       ' -v ',num2str(vmin), ...
                                       ' -i ',tmp_xdatafile, ...
                                       ' -j ',tmp_ydatafile, ...
                                       ' -o ',tmp_rqadatafile, ...
                                       ' -s']);
         errcode=22;
         % import RQA
         rqa_in = [];
         try
             fid = fopen(tmp_rqadatafile,'r'); % open RQA data file
             while 1
                 dataStr = fgetl(fid); % read RQA data
                 if ~ischar(dataStr), break, end % leave the loop if end of file
                 if isempty(findstr(dataStr, '#')) % neglect comments line (e.g. header)
                     rqa_in = str2num(dataStr); % import RQA measures into local variable rqa
                 end
             end
             fclose(fid); % close RQA data file

             RR = rqa_in(1);
             DET = rqa_in(2);
             LAM = rqa_in(4);
             Lmax = rqa_in(6);
             L = rqa_in(7);
             ENTR = rqa_in(8);
             Vmax = rqa_in(10);
             TT = rqa_in(11);
             t1 = rqa_in(13);
             t2 = rqa_in(14);
             RT = rqa_in(16);
             RTmax = rqa_in(15);
             RF = rqa_in(19);
             ENTW = rqa_in(17);
             warning off
             delete(tmp_rqadatafile);
             delete(tmp_xdatafile);
             delete(tmp_ydatafile);
             warning on
         catch 
             warning off
             delete(tmp_rqadatafile);
             delete(tmp_xdatafile);
             delete(tmp_ydatafile);
             warning on
         end

         if nogui~=2 & ishandle(h1), set(h1,'str',[num2str(i),'/',num2str(Nx-w)]); waitbar(i/(Nx-w)); drawnow, end



     % use builtin implementation
     else
     
         errcode=25;
         try
             if ~flag_pdist
                 if time_scale_flag
                   if length(x(i:i+w-1,:)) > 2000
                      X=crp_big(x(i:i+w-1,:),y(i:i+w-1,:),m,t,e,method,'nonorm','silent');
                   else
                      X=crp(x(i:i+w-1,:),y(i:i+w-1,:),m,t,e,method,'nonorm','silent');
                   end
                 else
                   X=crp2(x(i:i+w-1,:),y(i:i+w-1,:),m,t,e,method,'nonorm','silent');
                 end
           %  X=crp(x(i:i+w-1,:),y(i:i+w-1,:),m,t,e,varargin{i_char},'silent');
           else
               %%%%%%%%%%%%
               % alternative using pdist

               xcor = @(x,y) sqrt(sum((repmat(x,size(y,1),1)-y).^2,2));

               x_dist = pdist(x(i:i+w-1,:),xcor);
               X = squareform(x_dist);
           end
           
           %%%%%%%%%%%%
           
           warning off 
           if nogui~=2 & ishandle(h1), set(h1,'str',[num2str(i),'/',num2str(Nx-w)]); waitbar(i/(Nx-w)); drawnow, end

         %if 0
         catch
             error(lasterr)
             if nogui~=2 & ishandle(hw), close(hw), end
         end

         N=size(X);
         if theiler_window > 0
             X_theiler=double(triu(X,theiler_window))+double(tril(X,-theiler_window));
         else
             X_theiler=X;
         end
         errcode=26;

         % compute recurrence times of 1st and 2nd type
         if size(X_theiler,2) > 1000
             t1=[];t2=[];
             rps2=find(diff(double(X_theiler(:)))==1);
             rps=find(X_theiler(:));
             t1=diff(rps);
             t2=diff(rps2);
         else
             t1 = []; t2 = [];
             for i2=1:size(X_theiler,2)
                 if nogui~=2 & (Nx-w < 2) & ~rem(i2,25), waitbar(i2/size(X_theiler,2)), end
                 rps2=find(diff(double(X_theiler(:,i2)))==1);
                 rps=find(X_theiler(:,i2));
                 t1=[t1;diff(rps)];
                 t2=[t2;diff(rps2)];
             end
         end
         t1=mean(t1);
         t2=mean(t2);

         errcode=27;
         [a b]=dl(X_theiler);
         hist_l = [hist_l; b];
         N_hist_l(i) = length(b);

         warning off 
         errcode=271;
         b(b<lmin)=[];
         [c d dw]=tt(X_theiler);
         hist_v = [hist_v; d];
         hist_w = [hist_w; dw];         
         N_hist_v(i) = length(d);
         N_hist_w(i) = length(dw);
         
         warning off 
         errcode=272;
         d(find(d<vmin))=[];

         errcode=273;
         N_all = (N(1)*N(2));
         % reduce the number of possible states by the Theiler window;
         % because the Theiler window is applied symmetrically (on the LOI)
         % we only use N(1) and not N(2)
         if theiler_window >= 1
             N_all = N_all - N(1) - 2*((theiler_window-1)*N(1) - sum(1:(theiler_window-1)));
         end
         
         RR=sum(X_theiler(:))/N_all;
         %b(find(b>=max(N)-lmin))=[]; if isempty(b), b=0; end
         if isempty(b), b=0; end
         errcode=274;
         if sum(X_theiler(:)) > 0
           DET=sum(b)/sum(X_theiler(:));
         else
           DET=NaN;
         end
         errcode=275;
         L=mean(b);
         histL=hist(b(:),[1:min(N)]);
         ENTR=entropy(histL(:));
         errcode=276;
         if sum(X_theiler(:))>0
           LAM=sum(d)/sum(sum(X_theiler));
         else
           LAM=NaN;
         end
       
         % recurrence times
         RT = mean(dw);

         [dwh dwi] = hist(dw,[1:max(dw)]);
         if dwh
             [dws dwsi] = sort(dwh);
             RTp = dwi(dwsi(end)); % most probable recurrence time
         else
             RTp = 0;
         end
         RTmax = max(dw); % maximal recurrence time
         RF = 1/RTmax; % minimal recurrence frequency
         ENTW = entropy(dwh(:));

         errcode=277;
         TT=mean(d);
         b=[b;0]; Lmax=max(b);
         d=[d;0]; Vmax=max(d);
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %% new measures
         % RPDE
         errcode=278;
         bins = [.5:max(dw)+.5];
         
         dwh = histc(dw,bins); dwh = dwh(:);
         %dwh
         if max(dw) > 0 && numel(bins) > 3
            RTE = entropy(dwh) / log(max(dw));
         else
            RTE = NaN;
         end
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %% network measures
         % global clustering coefficient
         errcode=279;
         kv = sum(X_theiler,1); % degree of nodes
         Clust = mean(diag(double(X_theiler)*double(X_theiler)*double(X_theiler))' ./ (kv .* (kv-1)));
         % transitivity
         denom = sum(sum(double(X_theiler) * double(X_theiler)));
         Trans = trace(double(X_theiler)*double(X_theiler)*double(X_theiler))/denom;
         
     end % end plugin
  
     warning on

     errcode=28;
     Y(i,1)=RR; 
     Y(i,2)=DET;
     Y(i,3)=L;
     Y(i,4)=Lmax;
     Y(i,5)=ENTR;
     Y(i,6)=LAM;
     Y(i,7)=TT;
     Y(i,8)=Vmax;
     Y(i,9)=t1;
     Y(i,10)=t2;
     Y(i,11)=RTE;
     Y(i,12)=Clust;
     Y(i,13)=Trans;
     if undocumented
         Y(i,14)=RT;
         Y(i,15)=RTmax;
         Y(i,16)=RF;
         Y(i,17)=ENTW;
     else
         Y(i,15)=x_var;
         Y(i,16)=y_var;
         Y(i,17)=xy_var;
     end
 end % end window loop
 
% significance by bootstrap
 
% nboot = 100; 
% nD = round(mean(N_hist_l(N_hist_l>0)));
% nV = round(mean(N_hist_v(N_hist_v>0)));
% nW = round(mean(N_hist_w(N_hist_w>0)));
% for i = 1:nboot
%       onesample = ceil(length(hist_l)*rand(length(hist_l),1));
%       onesample = onesample(ceil(nD*rand(nD,1)));
%       tmp = sum(hist_l(onesample,:) >= lmin )/ sum(hist_l(onesample,:));
%       bootstatDET(i,:) = (tmp(:))';
%       tmp = mean(hist_l(onesample,:) >= lmin);
%       bootstatL(i,:) = (tmp(:))';
% 
%       onesample = ceil(length(hist_v)*rand(length(hist_v),1));
%       onesample = onesample(ceil(nV*rand(nV,1)));
%       tmp = sum(hist_v(onesample,:) >= vmin )/ sum(hist_v(onesample,:));
%       bootstatLAM(i,:) = (tmp(:))';
%       tmp = mean(hist_v(onesample,:) >= vmin);
%       bootstatTT(i,:) = (tmp(:))';
% 
%       onesample = ceil(length(hist_w)*rand(length(hist_w),1));
%       onesample = onesample(ceil(nW*rand(nW,1)));
%       tmp = mean(hist_w(onesample,:));
%       bootstatT2(i,:) = (tmp(:))';
% end

 
if ishandle(hw), waitbar(1); drawnow; close(hw), end % close waitbar

  if ~nogui
    h=findobj('tag','crqa_button_apply');
    set(h(1),'ToolTip','Starts the computation.','String','Apply','Callback','crqa compute')
    for j=1:length(obj); 
        h=findobj('Tag',obj{j},'Parent',h_fig(1)); 
        if ~isempty(h)
          set(h,'Enable','On')
        end
    end
  end

  if ~nogui, 
    h=findobj('Tag','crqa_Fig'); 
    if ~isempty(h), set(0,'CurrentFigure',h(1)); end
    setptr(gcf,'arrow'), 
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot
if ~nogui
    set(0,'showhidden','on')
    errcode=30;
    h=findobj('Tag','crqa_Fig'); if ~isempty(h), set(0,'CurrentFigure',h(1)); end
    tx={'RR';'DET';'L';'ENTR';'LAM';'TT';'RTE';'T_2';'Variance';'Covariance'};
    index=[1,2,3,5,6,7,9,11,15,16,17];
    tags={'crqa_axes_RR','crqa_axes_DET','crqa_axes_L','crqa_axes_ENTR','crqa_axes_LAM','crqa_axes_TT','crqa_axes_RTE','crqa_axes_T2','crqa_axes_Var','crqa_axes_CoVar'};
    h=findobj('Tag','crqa_axes_RR','Parent',gcf); h_axes.h(1)=h(1);
    h=findobj('Tag','crqa_axes_DET','Parent',gcf); h_axes.h(2)=h(1);
    h=findobj('Tag','crqa_axes_L','Parent',gcf); h_axes.h(3)=h(1);
    h=findobj('Tag','crqa_axes_ENTR','Parent',gcf); h_axes.h(4)=h(1);
    h=findobj('Tag','crqa_axes_LAM','Parent',gcf); h_axes.h(5)=h(1);
    h=findobj('Tag','crqa_axes_TT','Parent',gcf); h_axes.h(6)=h(1);
    h=findobj('Tag','crqa_axes_RTE','Parent',gcf); h_axes.h(7)=h(1);
    h=findobj('Tag','crqa_axes_T2','Parent',gcf); h_axes.h(8)=h(1);
    h=findobj('Tag','crqa_axes_Var','Parent',gcf); h_axes.h(9)=h(1);
    if ~all(x(:)==y(:)), h=findobj('Tag','crqa_axes_CoVar','Parent',gcf); h_axes.h(10)=h(1); end
   for i=1:9,
     set(gcf,'CurrentAxes',h_axes.h(i))
     if size(Y,1)==1
       cla
       text(0.5,0.5,sprintf('%6.4f',Y(index(i))),'FontWeight','bold','HorizontalAlign','Center')
%        if ~all(x(:)==y(:)) & i==9
%          cla 
%              text(0.5,0.5,sprintf('%6.5f, %6.5f',Y(index(i)),Y(index(i+1))),'FontWeight','bold','HorizontalAlign','Center')
%              set(gcf,'CurrentAxes',h_axes.h(i+1));cla
%          set(h_axes.h(i+1),'visible','off');
%        end
     else
       if i==9
         cla
             h2=stairs(xscale(1:wstep:length(Y)),Y(1:wstep:end,index(i)));
             set(h2,'color',props.line.Color)
         if ~all(x(:)==y(:))
           set(gca,'Color','none')
           hold on; h2=stairs(xscale(1:wstep:length(Y)),Y(1:wstep:end,index(10)),'r');
           set(gcf,'CurrentAxes',h_axes.h(i+1));cla
           set(h_axes.h(10),'visible','on');
           
           h3=stairs(xscale(1:wstep:length(Y)),Y(1:wstep:end,index(11))); set(h3,'color',[0 .4 0]);
           ylabel(tx(i+1));
           set(h_axes.h(i+1),'Tag',tags{i+1},'Units','Norm','Color','none','YAxisLocation','right','YColor',[0 .4 0])
         end
       else
         plot(xscale(1:wstep:length(Y)),Y(1:wstep:end,index(i)),'color',props.line.Color)
       end
     end
     set(gcf,'CurrentAxes',h_axes.h(i));
     ylabel(tx(i));
     set(gca,'Tag',tags{i},'color',props.axes.Color,'Units','Norm')
   end
   h=findobj('Tag','crqa_button_store');
   set(h(1),'UserData',Y,'Enable','On')
   xlim = get(h_axes.h(9),'xlim');
   if length(Y(:,1)) > 1, set(h_crqa_axes_Data,'xlim',xlim), end
   if ~all(x(:)==y(:)) & length(Y(:,1)) == 1, set(h_axes.h(10),'visible','off'); delete(get(h_axes.h(10),'children')), end
else
  if nargout, xout=Y(:,argout_index); end
end

if nargout, xout=Y(:,argout_index); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the end 

end




%if 0
catch
  if nogui~=2, if ishandle(hw), close(hw), end , end
  if nargout, xout = NaN; end
      
  z=whos;x=lasterr;y=lastwarn;in=varargin{1};
  print_error('crqa',z,x,y,in,method,action)
  try, if ~nogui
    h=findobj('tag','crqa_button_apply');
    set(h(1),'ToolTip','Starts the computation.','String','Apply','Callback','crqa compute')
    for j=1:length(obj); 
        h=findobj('Tag',obj{j},'Parent',h_fig(1)); 
        if ~isempty(h)
          set(h,'Enable','On')
        end
    end, end
  end
  set(0,'showhidden','on')
    h=findobj('Tag','crqa_Fig'); if ~isempty(h), setptr(h(1),'arrow'); end
  set(0,'showhidden','off')
end


try, set(0,props.root), end



try
  set(0,'ShowHidden','off')
end
