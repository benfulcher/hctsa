function X = taucrp(varargin)
%TAUCRP   Creates a close returns plot.
%    R=TAUCRP(X [,Y] [,param1,param2,...]) creates a cross 
%    recurrence plot/ recurrence plot R for a limited range of
%    past and future states, also known as close returns plot.
%
%    R=TAUCRP(X [,Y],M,T,E,W) uses the dimension M, delay T, the size 
%    of neighbourhood E and the range W of past and future time
%    steps.
%
%    If X and Y are multi-column vectors then they will be 
%    considered as phase space vectors (TAUCRP can be used
%    for real phase space vectors without embedding).
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
%      maxnorm     - Maximum norm, fixed recurrence rate.
%      fan         - Fixed amount of nearest neighbours.
%      distance    - Distance coded matrix (global CRP, Euclidean norm).
%
%    Normalization of the data series.
%      normalize   - Normalization of the data.
%      nonormalize - No normalization of the data.
%
%    Parameters not needed to be specified.
%
%
%    Examples: a = sin((1:1000) * 2 * pi/67);
%              w = 160;
%              X = taucrp(a,2,17,.2,w,'nonorm','euclidean');
%              imagesc(1:size(X,2),-w:w,X), colormap([1 1 1; 0 0 0])
%
%    See also CRP, CRP2, CRP_BIG, JRP, CRQA.
%
%    References: 
%    Marwan, N., Romano, M. C., Thiel, M., Kurths, J.: 
%    Recurrence Plots for the Analysis of Complex Systems, Physics 
%    Reports, 438(5-6), 2007.

% Copyright (c) 2008 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2008/01/25 12:47:26 $
% $Revision: 5.1 $
%
% $Log: taucrp.m,v $
% Revision 5.1  2008/01/25 12:47:26  marwan
% initial import
%


warning off
global errcode props nonorm

errcode=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% programme properties

init_properties
m_init = 1;
tau_init = 1;
eps_init = 0.1;
w_init = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check and read the input

error(nargchk(1,8,nargin));
if nargout>1, error('Too many output arguments'), end

check_meth={'ma','eu','mi','rr','fa','di'}; 	% maxnorm, euclidean, nrmnorm,  fan, distance
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

    if nonorm==1, 
	 x=(x0-repmat(mean(x0),Nx,1))./repmat(std(x0),Nx,1);
	 y=(y0-repmat(mean(y0),Ny,1))./repmat(std(y0),Ny,1);
    end

  if ~isempty(find(isnan(x))), for k=1:size(x,2),  x(find(isnan(x(:,k))),:)=[]; end, end
  if ~isempty(find(isnan(y))), for k=1:size(y,2),  y(find(isnan(y(:,k))),:)=[]; end, end
  if size(x,1) < t*(m0-1)+1 | size(y,1) < t*(m0-1)+1
     error(['Too less data',10,...
            'Either too much NaN or the number of columns in the vectors do not match.'])
  end

else
  error('No valid input given!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if m0>1
  x2=x(1:end-t*(m0-1),:);
  y2=y(1:end-t*(m0-1),:);
  for i=1:m0-1,
    x2(:,m*i+1:m*(i+1))=x(1+t*i:end-t*(m0-i-1),:);
    y2(:,m*i+1:m*(i+1))=y(1+t*i:end-t*(m0-i-1),:);
  end
  x=x2; y=y2; Nx=size(x,1); Ny=size(y,1);
  m=m0*m; clear x2 y2
end

x1=repmat(x,1,Ny);
for mi=1:m, x2(:,mi)=reshape(rot90(x1(:,0+mi:m:Ny*m+mi-m)),Nx*Ny,1); end
y1=repmat(y,Nx,1); x1=x2; clear x2


% create embedding vectors
NX = Nx - t*(m0-1); NY = Ny - t*(m0-1);

[idx add] = meshgrid(linspace(0,(m0-1)*t,m0), 1:Nx-(m0-1)*t);
idx = idx + add;
x1 = []; y1 = [];
for i = 1:m
    x1 = [x1 reshape(x(idx,i),size(idx,1),size(idx,2))];
    y1 = [y1 reshape(y(idx,i),size(idx,1),size(idx,2))];
end


% indices for which the recurrences should be computed
[i1 i2]  = meshgrid(1:NX,-w:w);
i2 = i2 + i1;

i_remove = i2 > length(y1) | i2 < 1;
i2(i_remove) = 1;

% compute distance matrix
D = (x1(i1,:) - y1(i2,:));

switch method
    case 1
        D2 = max(abs(D),[],2);
        X = double(D2 < e);
    case 2
        D2 = sqrt(sum(D.^2, 2));
        X = double(D2 < e);
    case 3
        D2 = sum(abs(D), 2);
        X = double(D2 < e);
    case 4
        D2 = max(abs(D),[],2);
        SS = sort(D2(:));
        idx = ceil(e * length(SS));
        e_ = SS(idx);
        X = double(D2 < e_);
    case 5
        D2 = sqrt(sum(D.^2, 2));
        D3 = reshape(D2,size(i1,1),size(i1,2));
        [SS, JJ] = sort(D3,1); JJ = JJ';
        mine = round((2*w+1)*e);
        X1(NX*(2*w+1)) = 0; 
        X1(JJ(:,1:mine)+repmat([0:(2*w+1):NX*(2*w+1)-1]',1,mine)) = 1;
        X = reshape(X1,(2*w+1),NX); 
    case 6
        D2 = sqrt(sum(D.^2, 2));
        X = double(D2);
      
    end

clear X1 SS JJ s px D D_ D2 D3
            
X = reshape(X,size(i1,1),size(i1,2)); X(i_remove) = 0;
%imagesc(X)
