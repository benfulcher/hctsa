function CPRout = phasesynchro(varargin);
% PS   Indicator of phase synchronisation by means of recurrences.
%    CPR=PHASESYNCHRO(X,Y [,param1,param2,...]) calculates the 
%    index of phase synchronisation based on recurrences.
%
%    CPR=PHASESYNCHRO(X,Y,M,T,E,W) uses the dimension M, delay T, 
%    the size of neighbourhood E and the range W of past and future 
%    time steps.
%
%    If X and Y are multi-column vectors then they will be 
%    considered as phase space vectors (TAUCRP can be used
%    for real phase space vectors without embedding).
%
%    The call of PHASESYNCHRO without output arguments plots the
%    tau-recurrence rate and the CPR value in the current figure.
%    
%    Parameters: dimension M, delay T, the size of neighbourhood 
%    E and the range W are the first four numbers after the data 
%    series; further parameters can be used to switch between 
%    various methods of finding the neighbours of the phasespace 
%    trajectory and to suppress the normalization of the data.
%
%    Methods of finding the neighbours/ of plot.
%      maxnorm     - Maximum norm.
%      euclidean   - Euclidean norm.
%      minnorm     - Minimum norm.
%      rr          - Maximum norm, fixed recurrence rate.
%      fan         - Fixed amount of nearest neighbours.
%
%    Normalization of the data series.
%      normalize   - Normalization of the data.
%      nonormalize - No normalization of the data.
%
%    Suppressing the plot of the results.
%      silent      - Suppresses the plot.
%
%    Parameters not needed to be specified.
%
%
%    Examples: a = sin((1:1000) * 2 * pi/67);
%              b = sin((1:1000) * 2 * pi/67) + randn(1,1000);
%              phasesynchro(a,b,2,17,'nonorm','euclidean');
%
%    See also CRP, CRP2, CRP_BIG, JRP, CRQA.
%
%    References: 
%    Marwan, N., Romano, M. C., Thiel, M., Kurths, J.: 
%    Recurrence Plots for the Analysis of Complex Systems, Physics 
%    Reports, 438(5-6), 2007.
%
%    Romano, M. C., Thiel, M., Kurths, J., Kiss, I. Z., Hudson, J.: 
%    Detection of synchronization for non-phase-coherent and 
%    non-stationary data, Europhysics Letters, 71(3), 2005.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2010/06/30 12:03:10 $
% $Revision: 5.6 $
%
% $Log: phasesynchro.m,v $
% Revision 5.6  2010/06/30 12:03:10  marwan
% Help text modified
%
% Revision 5.5  2009/06/02 14:07:39  marwan
% normalisation issue solved
%
% Revision 5.4  2009/03/24 08:34:24  marwan
% copyright address changed
%
% Revision 5.3  2008/07/01 11:36:05  marwan
% bug in embedding dimension fixed
%
% Revision 5.2  2008/04/29 14:49:30  marwan
% window size bug
%
% Revision 5.1  2008/01/25 12:47:25  marwan
% initial import
%
%
%
%    Examples: a = sin((1:1000) * 2 * pi/67);
%              b = sin((1:1000) * 2 * pi/67) + rand(1,1000);
%              X = phasesynchro(a,b,2,17,'nonorm','euclidean');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% programme properties

global errcode props

init_properties
m_init = 1;
tau_init = 1;
eps_init = 0.1;
w_init = 100;
method_str={'Maximum Norm','Euclidean Norm','Minimum Norm','Maximum Norm, fixed RR','FAN'};
norm_str = {'nor','non'};
lmin=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check the input

error(nargchk(1,9,nargin));
if nargout>1, error('Too many output arguments'), end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read the input


check_meth={'ma','eu','mi','rr','fa'}; 	% maxnorm, euclidean, nrmnorm,  fan
check_norm={'non','nor'};				% nonormalize, normalize
check_gui={'gui','nog','sil'};				% gui, nogui, silent

if isnumeric(varargin{1}) 		% read commandline input
   varargin{9}=[];
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

   % check the text input parameters for method, gui and normalization
   temp_meth=0;
   temp_norm=0;
   temp_gui=0;
   if ~isempty(i_char)
      for i=1:length(i_char), 
         varargin{i_char(i)}(4)='0';
         temp_meth=temp_meth+strcmpi(varargin{i_char(i)}(1:2),check_meth'); 
         temp_norm=temp_norm+strcmpi(varargin{i_char(i)}(1:3),check_norm'); 
         temp_gui=temp_gui+strcmpi(varargin{i_char(i)}(1:3),check_gui'); 
      end
      method=min(find(temp_meth));
      nonorm=min(find(temp_norm))-1;
      nogui=min(find(temp_gui))-1;
      if isempty(method), method=1; end
      if isempty(nonorm), nonorm=1; end
      if isempty(nogui), nogui=0; end
      if method>length(check_meth), method=length(check_meth); end
      if nonorm>1, nonorm=1; end
      if nogui>2, nogui=2; end
   else
      method=1; nonorm=1; nogui=0;
   end
   if nogui==0 & nargout>0, nogui=1; end

   % get the parameters for creating RP
     if max(size(varargin{1}))<=3
        error('To less values in data X.')
     end
     x=double(varargin{1});
     if isempty(varargin{2}) | ~isnumeric(varargin{2}), y=x; else
     y=double(varargin{2}); end

     if (isnumeric(varargin{2}) & max(size(varargin{2}))==1) | ~isnumeric(varargin{2})
       y=x;
       if ~isempty(varargin{i_double(2)}), m0=varargin{i_double(2)}(1); else m0=m_init; end
       if ~isempty(varargin{i_double(3)}), t=varargin{i_double(3)}(1); else t=tau_init; end
       if ~isempty(varargin{i_double(4)}), e=varargin{i_double(4)}(1); else e=eps_init; end
       if ~isempty(varargin{i_double(5)}), w=varargin{i_double(5)}(1); else w=w_init; end
     else
       if ~isempty(varargin{i_double(3)}), m0=varargin{i_double(3)}(1); else m0=m_init; end
       if ~isempty(varargin{i_double(4)}), t=varargin{i_double(4)}(1); else t=tau_init; end
       if ~isempty(varargin{i_double(5)}), e=varargin{i_double(5)}(1); else e=eps_init; end
       if ~isempty(varargin{i_double(6)}), w=varargin{i_double(6)}(1); else w=w_init; end
     end
     t=round(t); m0=round(m0); mflag=method;
     if e<0, e=1; disp('Warning: The threshold size E cannot be negative and is now set to 1.'), end
     if t<1, t=1; disp('Warning: The delay T cannot be smaller than one and is now set to 1.'), end
     if m0 < 1, m0 = 1; end
     if t < 1, t = 1; end
     if size(x,1)==1, x=x'; end, if size(y,1)==1, y=y'; end 
     m=max([size(x,2) size(y,2)]);
     if w > size(x,1); w = size(x,1); end
     if w > size(y,1); w = size(y,1); end

     if method==8 & (m*m0) > 1, 
       m0=1; 
       error(['The neighbourhood criterion ''Oder matrix''',10,'is not implemented - use crp or crp_big instead.'])
     end
     if method==9 & (m*m0) == 1, 
       m0=2; 
         disp(['Warning: For order patterns recurrence plots the dimension must',10,...
              'be larger than one. ',...
              'Embedding dimension is set to ',num2str(m0),'.'])
     end
     action='init';

  if ~isempty(find(isnan(x)))
     disp('NaN detected (in first variable) - will be cleared.')
     for k=1:size(x,2),  x(find(isnan(x(:,k))),:)=[]; end
  end
  if ~isempty(find(isnan(y)))
     disp('NaN detected (in second variable) - will be cleared.')
     for k=1:size(y,2),  y(find(isnan(y(:,k))),:)=[]; end
  end
  if size(x,1) < t*(m-1)+1 | size(y,1) < t*(m-1)+1
     error(['Too less data',10,...
            'Either too much NaN or the number of columns in the vectors do not match.'])
  end

    Nx=size(x,1); Ny=size(y,1);
    NX=Nx-t*(m0-1);NY=Ny-t*(m0-1);
    x0=zeros(Nx,m);y0=zeros(Ny,m);
    x0(1:size(x,1),1:size(x,2))=x; 
    y0(1:size(y,1),1:size(y,2))=y; 


  if ~isempty(find(isnan(x))), for k=1:size(x,2),  x(find(isnan(x(:,k))),:)=[]; end, end
  if ~isempty(find(isnan(y))), for k=1:size(y,2),  y(find(isnan(y(:,k))),:)=[]; end, end
  if size(x,1) < t*(m0-1)+1 | size(y,1) < t*(m0-1)+1
     error(['Too less data',10,...
            'Either too much NaN or the number of columns in the vectors do not match.'])
  end

else
  error('No valid input given!')
end




% compute RR-tau
X = taucrp(x, m0, t, e, w, check_meth(method),check_norm(nonorm+1));
X_ = taucrp(y, m0, t, e, w, check_meth(method),check_norm(nonorm+1));


%
X = [zeros(2*w+1,1) X zeros(2*w+1,1)];
X_ = [zeros(2*w+1,1) X_ zeros(2*w+1,1)];
x_diff = diff(X,[],2); % mark start and end points of line structures with -1 and +1
x_diff_ = diff(X_,[],2); % mark start and end points of line structures with -1 and +1

[line_start col1] = find(x_diff' == 1); % index of start points
[line_end col2] = find(x_diff' == -1); % index of end points
[line_start_ col1_] = find(x_diff_' == 1); % index of start points
[line_end_ col2_] = find(x_diff_' == -1); % index of end points


l = line_end - line_start; % length of all lines
l_corr = l; l_corr(l < lmin) = 0; % only lines longer than lmin
tau = unique(col1); % lag at which the lines where found
l_ = line_end_ - line_start_; % length of all lines
l_corr_ = l_; l_corr_(l_ < lmin) = 0; % only lines longer than lmin
tau_ = unique(col1_); % lag at which the lines where found
 

RR = zeros(2*w+1,1); RR_ = zeros(2*w+1,1);
for i = tau'; % waitbar(i/max(tau))
    j = col1==i;
    N_recpoints = sum(l(j)); % number of recurrence points at a given lag
    N_pointsline = sum(l_corr(j)); % number of recurrence points forming lines
    DET(i) = N_pointsline/N_recpoints;
    L(i) = mean(l_corr(j));
    RR(i) = sum(l(j))/(NX-abs(w-i));
end
for i = tau_';
    j = col1_==i;
    N_recpoints = sum(l_(j)); % number of recurrence points at a given lag
    N_pointsline = sum(l_corr_(j)); % number of recurrence points forming lines
    DET_(i) = N_pointsline/N_recpoints;
    L_(i) = mean(l_corr_(j));
    RR_(i) = sum(l_(j))/(NX-abs(w-i));
end
RR(find(isnan(RR)))=0; RR_(find(isnan(RR)))=0;



% normalise for covariance
% (start at w, i.e., at lag = 0)
rr1 = normalize(RR(w+1:end));
rr2 = normalize(RR_(w+1:end));

rr1a = (RR(w+1:end));
rr2a = (RR_(w+1:end));



% ACF of RR-tau
% find auto-correlation time
a1 = xcov(rr1,'unbiased'); a1(1:w) = [];
i1 = find(a1 < 1/exp(1),1);
a2 = xcov(rr2,'unbiased'); a2(1:w) = [];
i2 = find(a2 < 1/exp(1),1);
i3 = max([i1 i2]);


% correlation coefficient
C = corrcoef(rr1(i3:end), rr2(i3:end));
CPR = C(1,2);

% plot
if nogui ~= 2
    plot(0:i3,rr1a(1:i3+1),':',0:i3,rr2a(1:i3+1),':'), hold on
    plot(i3:w,rr1a(i3+1:end),i3:w,rr2a(i3+1:end)), hold off
    title(['CPR: ',sprintf('%2.2f',CPR)], 'fontw','bold')
    xlabel('Lag'), ylabel('RR_{\tau}')
end
if nargout | nogui == 2
    CPRout = CPR;
end
