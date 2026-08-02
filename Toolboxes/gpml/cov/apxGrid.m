function [K,Mx,xe] = apxGrid(cov, xg, hyp, x, z, b)

% apxGrid - Covariance function approximation based on an inducing point grid.
%
% A grid covariance function k(x,z) is composed as a product
% k(x,z) = k1(x(i1),z(i1)) * .. * kp(x(ip),z(ip)) of p covariance functions
% operating on mutually disjoint components of the data x,z.
% The resulting covariance matrix is given by the Kronecker product
%   K = kron(kron(kron(...,K3),K2),K1) = K_p x .. x K_2 x K_1.
%
% The function is designed to be used with infGaussLik and infLaplace.
%
% The covariance function grid contains p grid factors xg = {x1,..,xp}
% so that for each of the factors xi a covariance function ki can be defined.
% A factor xi is of size (ni,di) so that the Kronecker grid g has size
% (n1,n2,..,np,D), where D=d1+d2+..+dp. Hence, the grid g contains N=n1*n2*..*np
% data points overall. Note that the factors xi do neither need to
% be sorted nor do they need to be 1d a priori. If a factor xi contains unevenly
% spaced values, we require di=1.
%
% A factor xi can contain:
%  a) a univariate axis xi of size (ni,1),
%  b) a multivariate axis xi of size (ni,di),
%     where the values need to be equispaced, or
%  c) a stationary equispaced subgrid composed of di univariate equispaced axes
%     {xi_1,xi_2,xi_di}, where each axis xij is of size (nij,1) so that
%     ni=ni_1*..*ni_di intended to be used with stationary covariance functions
%     When a stationary equispaced subgrid is specified, we silently assume
%     the covariance function to be stationary.
%
% For fast computations, we exploit two kinds of structure:
%  1) The Kronecker structure of the covariance matrix induced BY the p factors.
%     Hence, for p=1, there is nothing to gain here.
%  2) The Toeplitz or BTTB (Block-Toeplitz with Toeplitz Blocks) WITHIN a factor
%     if a grid factor xi is equispaced (or multivariate and equispaced).
%     Note that empirically, only for matrices of sizes above 500x500, the
%     FFT-based MVMs are faster than dense matrix operations.
%
% Some examples with sizes and domain:
%  - a single factor with a univariate axis, case a)
%     xg = { 5*rand(150,1) }                => N = 150, D = 1, dom = [0,5]
%     xg = { linspace(0,3,400)' }           => N = 400, D = 1, dom = [0,3]
%  - a single factor with a multivariate axis (equispaced is mandatory), case b)
%     xg = { [linspace(0,3,100)',...
%             linspace(1,4,100)'] }         => N = 100, D = 2, dom = [0,3]x[1,4]
%  - a single factor with a univariate equispaced&stationary subgrid, case c)
%     where we assume a stationary covariance and exploit Toeplitz structure
%     xg = { {linspace(0,10,175)'} }        => N = 175, D = 1, dom = [0,10]
%  - a single factor with a bivariate equispaced&stationary subgrid, case c)
%     where we assume a stationary covariance and exploit BTTB structure
%     xg = { {linspace(0,3,20)',...
%             linspace(1,7,30)'} }          => N = 600, D = 2, dom = [0,3]x[1,7]
% - a two-factor grid of 2 univariate axes, case a)
%     xg = { 2*rand(40,1), 5*rand(20,1) }   => N = 800, D = 2, dom = [0,2]x[0,5]
% - a two-factor grid of 2 univariate equispaced and stationary axes, case c)
%     where we assume two stationary covariances and exploit Toeplitz structure
%     in both factors
%     xg = { {linspace(0,2,25)'}, ...
%            {linspace(1,3,25)'} }          => N = 625, D = 2, dom = [0,2]x[1,3]
% - a four-factor grid with a Toeplitz factor, two ordinary Kronecker
%   factors and a 2d BTTB factor
%     xg = { {linspace(0,1,50)'}, ...       => N = 4e6, D = 5, dom = [0,1]^5
%            rand(20,1), ...
%            linspace(0,1,40)', ... 
%            {linspace(0,1,10)',linspace(0,1,10)'} }
%
% The apxGrid function can be used to expand the (nested) cell array xg into a
% multivariate grid xe of size (N,D) via:
%     [xe,nx,Dx]  = apxGrid('expand',xg);                             => mode 1)
% The operation can be reverted (if no subgrids are used) by:
%     xg = apxGrid('factor',{xe,ng,Dg});                              => mode 2)
%
% Given scattered data x of size (n,D), we can create a grid xg covering
% the support of x using:
%     xg = apxGrid('create',x,eq,k);                                  => mode 3)
% The flag eq (default value 1) can be used to enforce an equispaced
% grid. The integer k scalar or vector, indicates the number of grid points per
% dimension. If k is a real number from (0,1], then the number of grid points
% equals k*numel(unique(x(:,1))).
% We require at least two different components per dimension.
%
% The variables v={x,z} can either be a) grid indices or b) data points.
% a) The variable v has size (nv,1) and contains integers from [1,N]. Then
%    the datapoints are obtained as g2 = reshape(g,N,D); v = g2(v,:).
% b) The variable v has size (nv,D) and directly represents the data points.
% The mechanism works for x and z separately.
%
% Given a grid xg, the grid size can be obtained by:
%     [ng,Dg] = apxGrid('size',xg);                                   => mode 4)
%
% An arbitrary data point x -- be it an index vector of size (n,1) or a data
% point of size (n,D) -- is converted into a regular data point xx of 
% size (n,D) by:
%     [xx,ng,Dg] = apxGrid('idx2dat',xg,x);                           => mode 5)
% If x is already of size (n,D), xx will simply equal x.
%
% Given a grid xg and given arbitrary data points x, the interpolation
% matrix Mx can directly be computed without computing cross-covariances via:
%     [Mx,dMx] = apxGrid('interp',xg,x,deg);                          => mode 6)
% For equispaced grids, deg can be used to set the degree of the interpolation
% polynomial in all p axes. Here deg=0 means nearest neighbor, deg=1 means
% linear interpolation, and deg=3 uses a cubic.
% The cell array dMx contains the derivatives d Mx / d xi with respect to
% the i=1..p grid components.
%
% Given a nested grid xg, we can compute a flattened grid xf by:
%     xf = apxGrid('flatten',xg);                                     => mode 7)
% without nesting and containing p axis elements.
%
% Given a covariance K, and a grid xg, and an integer embedding factor s,
% we can compute the Kronecker eigen decomposition such that
% K = V*diag(e)*V', where e = kron(ee,..).
%     [V,ee,e] = apxGrid('eigkron',K,xg,s);                           => mode 8)
% Note that this holds only approximately for Toeplitz/BTTB due to the
% circulant embedding. For details about s, see mode 9).
%
% Given a covariance function k where the call k(1:n) returns a vector of
% length, an interger length n and an integer embedding factor s, we
% construct a circulant embedding c (nx1).
%     [c,xg] = apxGrid('circ',k,n,s);                                 => mode 9)
% Here s is allowed to have the following values:
%     s=0  No embedding c = k.
%     s>0  Whittle embedding [1] as described by Guinness & Fuentes [2] in
%          equation (5) with N = |s|.
%     [1] Whittle, On stationary processes in the plane, Biometrika, 1954.
%     [2] Guinness & Fuentes, Circulant embedding of approximate covariances for
%         inference from Gaussian data on large lattices, 2014.
%
% Given a grid covariance K, a grid xg, an interpolation matrix Mx and 
% two sets of column vectors a and b, we can compute the directional derivative
% d trace(a'*Mx*K*Mx'*b) / d hyp:
%     dhyp = apxGrid('dirder',K,xg,Mx,a,b);                          => mode 10)
%
% Multiplication with a matrix/operator A along dimension dim of the tensor B.
% Note that tmul(A,B,1) = A*B if B is a matrix.
%     C = apxGrid('tmul',A,B,dim);                                   => mode 11)
%
% Test grid for being equispaced along an axis i
%     eq = apxGrid('equi',xg,i);                                     => mode 12)
%
% Return a descriptive string about the nature of Kg and Mx.
%     s = apxGrid('info',Kg,Mx,xg,deg);                              => mode 13)
%
% PerformLanczos with at most d MVMs with B.
% [Q,T] = apxGrid('lanczos_arpack',B,v,d)                            => mode 14)
%
% The hyperparameters are:
% hyp = [ hyp_1
%         hyp_2
%          ..
%         hyp_p ],
%
% Copyright (c) by Hannes Nickisch, Kun Dong, Insu Han 
%                                                  and Andrew Wilson 2018-03-26.
%
% See also covFunctions.m, cov/apx.m, inf/infLaplace.m, inf/infGaussLik.m.

if nargin<2, error('Not enough parameters provided.'), end
dense_max = 0; % 500 if larger grid we use FFT-algebra rather than dense algebra

% mode  1) expand axes xg representation into full grid x
if     strcmp(cov,'expand')          % call: [xe,nx,Dx]  = apxGrid('expand',xg);
  [K,Mx,xe] = expandgrid(xg); return

% mode  2) factor full x grid into axes representation xg
elseif strcmp(cov,'factor')           % call: xg = apxGrid('factor',{xe,ng,Dg});
  K = factorgrid(xg{:}); return

% mode  3) create axes representation xg from scattered data
elseif strcmp(cov,'create')               % call: xg = apxGrid('create',x,eq,k);
  if nargin<3, eq = 1; else eq = hyp; end                   % set default values
  if nargin<4, k = 1; else k = x; end, x = xg;                % set input params
  K = creategrid(x,eq,k); return

% mode  4) convert possible index vector into data space
elseif strcmp(cov,'size')                  % call: [ng,Dg] = apxGrid('size',xg);
  [ng,Dg] = sizegrid(xg);
  K = ng; Mx = Dg; return

% mode  5) convert possible index vector into data space
elseif strcmp(cov,'idx2dat')       % call: [xx,ng,Dg] = apxGrid('idx2dat',xg,x);
  [ng,Dg] = sizegrid(xg); N = prod(ng); 
  if isidx(hyp,N), xe = expandgrid(xg); K = xe(hyp,:); else K = hyp; end
  Mx = ng; xe = Dg; return

% mode  6) compute interpolation matrix
elseif strcmp(cov,'interp')        % call [Mx,dMx] = apxGrid('interp',xg,x,deg);
  if nargout>1, [K,Mx]=interpgrid(xg,hyp,x); else K=interpgrid(xg,hyp,x); end
  return

% mode  7) provide flattened interpolation grid without nesting
elseif strcmp(cov,'flatten')                  % call xf = apxGrid('flatten',xg);
  K = flattengrid(xg); return

% mode  8) compute eigen-decomposition of Kronecker matrix
elseif strcmp(cov,'eigkron')        % call [V,ee,e] = apxGrid('eigkron',K,xg,s);
  [K,Mx,xe] = eigkron(xg,hyp,x); return

% mode  9) compute a circulant embedding
elseif strcmp(cov,'circ')                 % call [c,xg] = apxGrid('circ',k,n,s);
  [K,Mx] = circ(xg,hyp,x); return

% mode 10) compute directional derivatives
elseif strcmp(cov,'dirder')         % call dhyp = apxGrid('dirder',K,xg,Mx,a,b);
  K = dirder(xg,hyp,x,z,b); return

% mode 11) matrix multiplication along a tensor
elseif strcmp(cov,'tmul')                    % call C = apxGrid('tmul',A,B,dim);
  K = tmul(xg,hyp,x); return

% mode 12) test grid for being equispaced
elseif strcmp(cov,'equi')                      % call eq = apxGrid('equi',xg,i);
  K = equi(xg,hyp); return 

% mode 13) report a descriptive status
elseif strcmp(cov,'info')               % call s = apxGrid('info',Kg,Mx,xg,deg);
  K = info(xg,hyp,x,z); return

% mode 14) perform Lanczos
elseif strcmp(cov,'lanczos_arpack')% call [Q,T]=apxGrid('lanczos_arpack',B,v,d);
  [K,Mx] = lanczos_arpack(xg,hyp,x); return
end

% mode  0) regular covariance function computations
p = numel(xg); [ng,Dg] = sizegrid(xg);             % number of Kronecker factors
if numel(cov)~=p, error('We require p factors.'), end
for ii = 1:p                                 % iterate over covariance functions
  f = cov(ii); if iscell(f{:}), f = f{:}; end   % expand cell array if necessary
  D = Dg(ii); j(ii) = cellstr(num2str(eval(feval(f{:}))));  % collect nbr hypers
end

if nargin<4                                        % report number of parameters
  K = char(j(1)); for ii=2:length(cov), K = [K, '+', char(j(ii))]; end, return
end
if nargin<5, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

v = [];               % v vector indicates to which covariance parameters belong
for ii = 1:p, v = [v repmat(ii, 1, eval(char(j(ii))))]; end

N = prod(ng); n = size(x,1); D = sum(Dg);     % expanded grid and data dimension
ix = isidx(x,N);               % determine whether x is an index or a data array
if ~ix && size(x,2)~=D, error('Grid and data dimension are different.'), end
if nargout>1 || ~(dg||xeqz||ix)                         % off-grid interpolation
  if nargin>5, Mx = interpgrid(xg,x,b); else Mx = interpgrid(xg,x); end
end
if dg               % evaluate as full dense vector for diagonal covariance case
  K = 1;                       % xg is not assumed to form a grid for z = 'diag'
  for ii = 1:length(cov)                       % iteration over factor functions
    f = cov(ii); if iscell(f{:}), f = f{:}; end % expand cell array if necessary
    d = sum(Dg(1:ii-1))+(1:Dg(ii));                     % dimensions of interest
    if ix, xii = xg{ii}; else xii = x(:,d); end    % switch Kronecker/plain prod
    Kj = feval(f{:}, hyp(v==ii), xii, z);               % plain Kronecker factor
    if ix, K = kron(K,Kj); else K = K.*Kj; end     % switch Kronecker/plain prod
  end
  if ix, K = K(x); end, return
end

if isidx(z,N), iz = z; z = apxGrid('expand',xg); z = z(iz,:); end     % expand z
K = cell(p,1); dK = cell(p,1); sz = [1,1];    % cov Kronecker factors total size
for ii = 1:p                                   % iteration over factor functions
  f = cov(ii); if iscell(f{:}), f = f{:}; end   % expand cell array if necessary
  d = sum(Dg(1:ii-1))+(1:Dg(ii));                       % dimensions of interest
  if isnumeric(z) && ~isempty(z)                                   % cross terms
    zd = z(:,d);
  elseif xeqz && iscell(xg{ii})                                  % Toeplitz/BTTB
    zd = [];
  else                                                        % symmetric matrix
    zd = z;
  end
  xii = xg{ii};                                                    % grid factor
  if xeqz && iscell(xg{ii})                                      % Toeplitz/BTTB
    di = numel(xii); ni = sizegrid(xii);
    wi = zeros(di,1); for j=1:di, wi(j) = xii{j}(end)-xii{j}(1); end     % width
    eqstat = true; for j=1:di, eqstat = eqstat & equi(xii,j); end
    if ~eqstat, error('Subgrid not equispaced.'), end          % stop if problem
    if prod(ni)>dense_max    % empirical thresh: FFT-based MVM with 1 rhs faster
      kii = @(n) feval(f{:},hyp(v==ii),(n-1)*diag(wi./(ni-1)), zeros(1,di));
      xc = cell(di,1);  % generic (Strang) circular embedding grid, kii = k(int)
      for j=1:di, n2 = floor(ni(j)-1/2)+1; xc{j} = [1:n2,n2-2*ni(j)+2:0]'; end
      ci = kii(apxGrid('expand',xc)); ci = reshape(ci,[2*ni(:)'-1,1]); dki = [];
      fi = real(fftn(ci));                  % precompute FFT for circular filter
      mvmKi = @(x) bttbmvmsymfft(fi,x);          % MVM with Toeplitz/BTTB matrix
      if di==1, s = 'toep'; else s = ['bttb',num2str(di)]; end
      ki = struct('descr',s, 'mvm',mvmKi, 'kii',kii, 'size',prod(ni)*[1,1]);
    else               % simply evaluate covariance matrix if prod(ni) too small
      [ki,dki] = feval(f{:},hyp(v==ii),expandgrid(xii));
    end
    K{ii} = ki; dK{ii} = dki; sz = sz.*K{ii}.size;
  else
    if iscell(xii), xii = expandgrid(xii); end
    [K{ii},dK{ii}] = feval(f{:}, hyp(v==ii), xii, zd);  % plain Kronecker factor
    sz = sz.*size(K{ii});
  end
end

if xeqz                                                    % create mvm and rest
  K = struct('mvm',@(a)kronmvm(K,a),'size',sz,'kronmvm',@kronmvm,...
                                        'kron',struct('factor',K,'dfactor',dK));
else                                                        % expand cross terms
  Ks = K; K = Ks{1}; for ii = 2:p, K = kron1(Ks{ii},K); end
  if ix, if numel(x)~=N || max(abs(x-(1:N)'))>0, K = K(x,:); end
  else   K = Mx*K;
  end
end
if nargout>2, xe = apxGrid('expand',xg); end

% Perform a matrix vector multiplication b = A*x with a matrix A being a
% Kronecker product given by A = kron( kron(...,As{2}), As{1} ).
function b = kronmvm(As,x,transp)
if nargin>2 && ~isempty(transp) && transp   % transposition by transposing parts
  for i=1:numel(As)
    if isnumeric(As{i})
      As{i} = As{i}';
    else
      As{i}.mvm = As{i}.mvmt;
      As{i}.size = [As{i}.size(2),As{i}.size(1)];
    end
  end
end
m = zeros(numel(As),1); n = zeros(numel(As),1);                  % extract sizes
for i=1:numel(n)
  if isnumeric(As{i})
    [m(i),n(i)] = size(As{i});
  else
    m(i) = As{i}.size(1); n(i) = As{i}.size(2);
  end
end
d = size(x,2);
b = x;
for i=1:numel(n)                              % apply As{i} to the 2nd dimension
  sa = [prod(m(1:i-1)), n(i), prod(n(i+1:end))*d];                        % size
  a = reshape(permute(reshape(full(b),sa),[2,1,3]),n(i),[]);
  if isnumeric(As{i}), b = As{i}*a; else b = As{i}.mvm(a); end    % do batch MVM
  b = permute(reshape(b,m(i),sa(1),sa(3)),[2,1,3]);
end
b = reshape(b,prod(m),d);                        % bring result in correct shape

% Perform MVM b = T*a with a of size (n,m) with a BTTB (Block-Toeplitz with
% Toeplitz-blocks) matrix T of size (n,n) by pointwise multiplication with the
% Fourier-transformed filter f. All variables are assumed real valued.
% Needs O(3*m*n*log(n)) time and O(n*m) space.
function b = bttbmvmsymfft(f,a)
  ng = (size(f)+1)/2; p = numel(ng); Ng = prod(ng);              % extract sizes
  if p==2 && ng(2)==1, p = 1; ng = ng(1); end           % detect 1d and reduce p
  m = numel(a)/Ng; b = reshape(a,[ng,m]);
  for i=1:p, b = fft(b,2*ng(i)-1,i); end           % emulate fftn with new shape
  b = bsxfun(@times,f,b);                             % pointwise multiplication
  for i=1:p, b = ifft(b,[],i); end                                % emulate ifft
  for i=1:p                          % only keep the relevant part of the result
    b = reshape(b,prod(ng(1:i-1)),2*ng(i)-1,prod(2*ng(i+1:p)-1)*m);
    b = b(:,1:ng(i),:);
  end
  b = real(reshape(b,[],m));

% perform kron along first dimension only
% the code is equivalent to the following loop
%   z = zeros(size(x,1)*size(y,1),size(x,2));
%   for i=1:size(z,2), z(:,i) = kron(x(:,i),y(:,i)); end
function z = kron1(x,y)
  nx = size(x,1); ny = size(y,1);
  z = repmat(reshape(x,1,nx,[]),[ny,1,1]).*repmat(reshape(y,ny,1,[]),[1,nx,1]);
  z = reshape(z,nx*ny,[]);

function r = isidx(i,N)     % check whether i represents an integer index vector
  r = false;
  if numel(i)>0 && ~strcmp(i,'diag') && size(i,2)==1 && ndims(i)==2
    if max(abs(i-floor(i)))<1e-13
      if 0<min(i) && max(i)<=N, r = true; end
    end
  end
  
% mode 1 (expand)
function [x,ng,Dg] = expandgrid(xg)                    % expand a Kronecker grid
  [ng,Dg] = sizegrid(xg);                                        % original size
  if ~iscell(xg), x = xg; return, end                       % catch trivial case
  xg = flattengrid(xg);                                      % remove nestedness
  p = numel(xg); x = xg{1};                                 % expanded grid data
  ngf = zeros(p,1); Dgf = zeros(p,1); [ngf(1),Dgf(1)] = size(xg{1});
  for i=2:p
    szx = size(x); [ngf(i),Dgf(i)] = size(xg{i});
    xold = repmat(reshape(x,[],1,szx(end)),[1,ngf(i),1]);
    xnew = repmat(reshape(xg{i},[1,ngf(i),Dgf(i)]),[size(xold,1),1,1]);
    x = reshape(cat(3,xold,xnew),[szx(1:end-1),ngf(i),szx(end)+Dgf(i)]);
  end
  x = reshape(x,[],size(x,ndims(x)));

% mode 2 (factor)
function xg = factorgrid(x,ng,Dg)                      % factor a Kronecker grid
  p = numel(ng); xg = cell(p,1);         % extract individual grid components xg
  for i=1:p
    x = reshape(x,[prod(ng(1:i-1)), ng(i), prod(ng(i+1:end)), sum(Dg)]);
    xg{i} = reshape(x(1,:,1,sum(Dg(1:i-1))+(1:Dg(i))), ng(i), Dg(i));
  end

% mode 3 (create)
function xg = creategrid(x,eq,k)
  if nargin<2, eq = 1; end                                  % set default values
  if nargin<3, k = 1; end                                     % set input params
  p = size(x,2); xg = cell(p,1);                               % allocate result
  if numel(k)>0, k = ones(p,1).*k(:); end              % enforce vector-valued k
  for j=1:p                                            % iterate over dimensions
    u = sort(unique(x(:,j))); if numel(u)<2, error('Two few unique points.'),end
    if isempty(k)                              % determine number of grid points
      if eq
        ngj = ceil( (u(end)-u(1))/min(abs(diff(u))) );     % use minimum spacing
      else
        ngj = numel(u);
      end
    elseif 0<=k(j) && k(j)<=1
      ngj = ceil(k(j)*numel(u));
    else
      ngj = k(j);
    end
    du = (u(end)-u(1))/ngj; bu = [u(1)-5*du, u(end)+5*du];
    if eq                                                      % equispaced grid
      xg{j} = linspace(bu(1),bu(2),max(ngj,5))';        % at least 5 grid points
    else                                                   % non-equispaced grid
      [idx,xgj] = kmeans(u,min(numel(u),ngj-2)); xgj = sort(xgj(:))';  % cluster
      nb = ngj-numel(xgj); nb1 = floor(nb/2); nb2 = nb - nb1; % size of boundary
      xg1 = linspace(bu(1),xgj(1),nb1+1); xg2 = linspace(xgj(end),bu(2),nb2+1);
      xg{j} = [xg1(1:nb1), xgj, xg2(1+(1:nb2))]';
    end
  end

% mode 4 (size)
function [ng,Dg] = sizegrid(xg)          % report the size of the p grid factors
  if ~iscell(xg), [ng,Dg] = size(xg); return, end           % catch trivial case
  p = numel(xg); ng = zeros(p,1); Dg = zeros(p,1);      % number of grid factors
  for i=1:p                                          % iterate over grid factors
    x = xg{i};
    if iscell(x)                                 % stationary and equispace grid
      for j=1:numel(x)
        if j==1
          [ng(i),Dg(i)] = size(x{j});
        else
          ng(i) = ng(i)*size(x{j},1);
          Dg(i) = Dg(i)+size(x{j},2);
        end
      end
    else                                                        % arbitrary grid
      [ng(i),Dg(i)] = size(x);
    end
  end

% mode 6 (interp)
% deg, degree of equispaced interpolation polynomial, 0:nn, 1:lin, 3:cub
function [Mx,dMx] = interpgrid(xg,x,deg)
  xg = flattengrid(xg);                                      % remove nestedness
  p = numel(xg); Dg = zeros(p,1); ng = zeros(p,1); n = size(x,1);      % dims ..
  for i=1:p, [ng(i),Dg(i)] = size(xg{i}); end, N = prod(ng);       %.. and sizes
  ix = isidx(x,N);             % determine whether x is an index or a data array
  if ix
    Mx = sparse(1:n,x,1,n,N); if nargout>1, dMx = repmat({sparse(n,N)},p,1); end
  else
    if nargin<3, deg = 3; end, deg = deg(:).*ones(p,1);       % cubic is default
    s = 1;                                                      % initial stride
    for i=1:p                                 % iterate over Toeplitz components
      d = sum(Dg(1:i-1))+(1:Dg(i));                     % dimensions of interest
      xt = xg{i}; it = find(abs(xt(2,:)-xt(1,:)));        % grid nonzero inc idx
      if equi(xg,i)                         % compute interpolation coefficients
        if nargout>1                                       % equispaced grid pts
          [Ji,Ci,dCi] = eqinterp(xt(:,it),x(:,d(it)),deg(i));
        else
          [Ji,Ci] = eqinterp(xt(:,it),x(:,d(it)),deg(i));
        end
      else
        if nargout>1   % non-equispaced grid pts, lin interp, inv dist weighting
          [Ji,Ci,dCi] = neqinterp(xt(:,it),x(:,d(it)));
        else
          [Ji,Ci] = neqinterp(xt(:,it),x(:,d(it)));
        end
      end
      nc = size(Ci,2);    % number of interpolation coefficients along dimension
      if i==1
        C = Ci; J = ones(n,1);
        if nargout>1, dC = repmat({Ci},p,1); dC{1} = dCi; end
      else
        C = repmat(C,[1,1,nc]) .* repmat(reshape(Ci,n,1,nc),[1,size(C,2),1]);
        C = reshape(C,n,[]);
        if nargout>1
          for j=1:p
            if i==j, dCij = dCi; else dCij = Ci; end
            dC{j} = repmat(dC{j},[1,1,nc]) .* ...
                    repmat(reshape(dCij,n,1,nc),[1,size(dC{j},2),1]);
            dC{j} = reshape(dC{j},n,[]);
          end
        end
      end
      J = repmat(J(:),[1,nc]) + s*repmat(Ji-1,[size(C,2)/nc,1]);  % blow 2nd idx
      s = s*ng(i);                                               % update stride
    end
    I = repmat((1:n)',[1,size(C,2)]); id = 0<J&J<=N;% first index and valid flag
    Mx = sparse(I(id),J(id),C(id),n,N);
    if nargout>1
      dMx = cell(p,1); for i=1:p,dMx{i} = sparse(I(id),J(id),dC{i}(id),n,N); end
    end
  end

% Compute interpolation coefficients C (nt,nc) and interpolation coefficient
% indices J (nt,nc) from a source grid s (ns,1) to a target array t (nt,1).
% The coefficient matrix C has rows summing up to 1.
function [J,C,dC] = eqinterp(s,t,d)
  gp = false;
  switch d
    case 0, k = @(x) -0.5<x & x<=0.5; it=-1:0; dk = @(x) 0*x;
    case 1, k = @(x) max(1-abs(x),0); it=-1:0; dk = @(x) -(abs(x)<=1).*sign(x);
    case 3, k = @kcub;                it=-2:1; dk = @dkcub;
    otherwise, ell = d/5; gp = true;                          % GP interpolation
      k = @(x) exp(-x.*x/(2*ell^2));  it=(0:d-1)-floor(d/2);
      dk = @(x) -k(x).*x/(ell^2);
  end
  ds = s(2)-s(1); ns = numel(s); nt = numel(t); nc = numel(it);
  if size(s,2)*size(t,2)~=1, error('Interpolation only possible for d==1.'), end
  if ns<nc, error('Interpolation only possible for ns>%d.',nc-1), end
  j = floor((t-s(1))/ds)+1;                % index of closest smaller grid point
  w = (t-s(1))/ds-j+1;   % relative distance to closest smaller grid point [0,1]
  j = j-it(nc); C = zeros(nt,nc); dC = zeros(nt,nc);
  for i=1:nc
    C(:,i) = k(w+it(nc+1-i)); if nargout>2, dC(:,i)=dk(w+it(nc+1-i))*(ns-1); end
  end
  if gp, kn = k(sqrt(sq_dist(1:d))); C = C/kn; dC = dC/kn; end
  v = 1; id = find(j<nc+it(1)); C(id,:) = 0; dC(id,:) = 0;  % fix lower boundary
  D = abs(repmat(s(1:nc)',numel(id),1)-repmat(t(id),[1,nc]));
  [junk,jid] = min(D,[],2);          % index of closest index in boundary region
  for i=1:numel(id), C(id(i),jid(i)) = 1; dC(id(i),jid(i)) = 0; end, j(id) = v;
  v = ns-nc+1; id = find(j>v); C(id,:) = 0; dC(id,:) = 0;   % fix upper boundary
  D = abs(repmat(s(ns-nc+1:ns)',numel(id),1)-repmat(t(id),[1,nc]));
  [junk,jid] = min(D,[],2);          % index of closest index in boundary region
  for i=1:numel(id), C(id(i),jid(i)) = 1; dC(id(i),jid(i)) = 0; end, j(id) = v;
  J = zeros(nt,nc); for i=1:nc, J(:,i) = j+i-1; end    % construct index array J

% Robert G. Keys, Cubic Convolution Interpolation for Digital Image Processing,
% IEEE ASSP, 29:6, December 1981, p. 1153-1160.
function y = kcub(x)
  y = zeros(size(x)); x = abs(x);
  q = x<=1;          % Coefficients:  1.5, -2.5,  0, 1
  y(q) =            (( 1.5 * x(q) - 2.5) .* x(q)    ) .* x(q) + 1;
  q = 1<x & x<=2;    % Coefficients: -0.5,  2.5, -4, 2
  y(q) =            ((-0.5 * x(q) + 2.5) .* x(q) - 4) .* x(q) + 2;
function y = dkcub(x)
  y = sign(x); x = abs(x);
  q = x<=1;          % Coefficients:  1.5, -2.5,  0, 1
  y(q) = y(q) .*  ( 4.5 * x(q) - 5.0) .* x(q);
  q = 1<x & x<=2;    % Coefficients: -0.5,  2.5, -4, 2
  y(q) = y(q) .* ((-1.5 * x(q) + 5.0) .* x(q) - 4.0);
  y(x>2) = 0;

% Perform piecewise linear interpolation using inverse distance weighting.
% s (ns,1) source nodes, need neither be sorted nor equispaced
% t (nt,1) target nodes
% M (nt,ns) interpolation matrix, M = sparse((1:N)'*[1,1],J,C,nt,ns);
%
% z = M*y where y (ns,1) are source values and z (nt,1) are target values
function [J,C,dC] = neqinterp(s,t)
  ns = size(s,1); nc = 2;                                       % get dimensions
  if size(s,2)*size(t,2)~=1, error('Interpolation only possible for d==1.'), end
  if ns<nc, error('Interpolation only possible for ns>=nc.'), end
  [s,ord] = sort(s); ds = diff(s); 
  if min(ds)<1e-10, error('Some source points are equal.'), end
  [junk,ii] = histc(t(:),[-inf;s(2:end-1);inf]);
  d0 = t(:)-s(ii); d1 = s(ii+1)-t(:); d0n = d0<0; d1n = d1<0;
  d0(d0n) = 0; d1(d1n) = 0;                                % boundary conditions
  J = [ord(ii),ord(ii+1)]; C = [d1./(d1+d0),d0./(d1+d0)];
  nz = 1-(d1n|d0n); if nargout>2, dC = [-nz./(d1+d0),nz./(d1+d0)]; end

% mode 7 (flatten)
function xf = flattengrid(xg)               % convert nested grid into flat grid
  if ~iscell(xg), xf = xg; return, end                      % catch trivial case
  xf = cell(1,0);
  for i=1:numel(xg)
    x = xg{i}; if iscell(x), xf = [xf,x]; else xf = [xf,{x}]; end
  end

% mode 8 (eigkron)
% Eigendecomposition of a Kronecker matrix K with dense, Toeplitz or BTTB
% factors so that K = V*diag(e)*V', where e = kron(ee,..). Note that this holds
% only approximately for Toeplitz/BTTB due to the circulant embedding.
function [V,ee,e] = eigkron(K,xg,s)
isbttb = @(Ki) isstruct(Ki) && (strcmp (Ki.descr,'toep') ...
                             || strncmp(Ki.descr,'bttb',4));   % BTTB covariance
p = numel(K.kron); V = cell(p,1); ee = cell(p,1);    % sizes and allocate memory
for j=1:p                                   % compute eigenvalue diagonal matrix
  if isbttb(K.kron(j).factor)
    [xj,nj] = apxGrid('expand',xg{j});                    % extract subgrid size
    ej = fftn(circ(K.kron(j).factor.kii,nj,s));          % circ embedded cov mat
    V{j}.mvm  = @(v) Fmvm( v,nj);                      % V{j}' is Fourier matrix
    V{j}.mvmt = @(v) Fmvmt(v,nj); V{j}.size = numel(ej)*[1,1];
  else
    Kj = K.kron(j).factor; Kj = (Kj+Kj')/2;                   % enforce symmetry
    [V{j},ej] = eig(Kj);            % compute eigenvalues of non-Toeplitz matrix
    ej = diag(ej); V{j} = real(V{j});                                % cosmetics
  end
  ee{j} = max(real(ej),0);             % thresholding to ensure pd approximation
end
if nargout>2, e = 1; for j=1:p, e = kron(ee{j}(:),e); end, end

function a = Fmvmt(b,nj)     % fast Fourier transform transpose for multiple rhs
  Nj = prod(nj); sNj = sqrt(Nj);       % scaling factor to make FFTN orthonormal
  nr = numel(b)/Nj;                        % number of right-hand-side arguments
  b = reshape(b,[nj(:)',nr]);
  for i=1:numel(nj), b = fft(b,[],i); end                         % emulate fftn
  a = reshape(b,Nj,[])/sNj;                                  % perform rescaling

function b = Fmvm(a,nj)                % fast Fourier transform for multiple rhs
  Nj = prod(nj); sNj = sqrt(Nj);       % scaling factor to make FFTN orthonormal
  nr = numel(a)/Nj;                        % number of right-hand-side arguments
  b = a; b = reshape(b,[nj(:)',nr]);               % accumarray and target shape
  for i=1:numel(nj), b = ifft(b,[],i); end                       % emulate ifftn
  b = reshape(b,Nj,nr)*sNj;                                  % perform rescaling

% mode 9 (circ)
% Construct a circular embedding c(nx1) from a covariance function k.
%  - k is a function and the call k(1:n) returns a vector of length n
%  - s is the setting for the embedding with values s={0,1,2,..}.
%
% s=0  No embedding c = k.
% s>0  Whittle embedding [1] as described by Guinness & Fuentes [2] in
%      equation (5) with N = |s|.
%
% [1] Whittle, On stationary processes in the plane, Biometrika, 1954, 41(3/4).
% [2] Guinness & Fuentes, Circulant embedding of approximate covariances for
%     inference from Gaussian data on large lattices, 2014, preprint,
%     http://www4.stat.ncsu.edu/~guinness/circembed.html.
function [c,xg] = circ(k,n,s)
p = numel(n); n = n(:)';                             % dimensions and row vector
if nargin<3, s = 2; end                                          % default value
if s==0                                                    % no embedding at all
  xg = cell(p,1); for i=1:p, xg{i} = (1:n(i))'; end              % standard grid
  c = reshape(k(apxGrid('expand',xg)),[n,1]);
elseif s>0                                           % Whittle/Guinness aliasing
  xg = cell(p,1); for i=1:p, xg{i} = (1-s*n(i):s*n(i))'; end
  sz = [n; 2*s*ones(1,p)]; c = reshape(k(apxGrid('expand',xg)), sz(:)');
  for i=1:p, c = sum(c,2*i); end, c = squeeze(c);
end

% mode 10 (dirder)
% d trace(a'*Mx*Kg*Mx'*b) / d hyp
function dhyp = dirder(Kg,xg,Mx,a,b)
  a = reshape(a,size(Mx,1),[]);                       % turn into column vectors
  if nargin<5, b = a; end                        % default input if b is missing
  b = reshape(b,size(Mx,1),[]);                       % turn into column vectors
  p = numel(Kg.kron); na = size(a,2);     % number of Kronecker factors, vectors
  ng = [apxGrid('size',xg)',na];                                % grid dimension
  Mta  = Mx'*a; Mtb = Mx'*b; dhyp = [];              % dhyp(i) = trace(a'*dKi*b)
  for i=1:p
    sz = [prod(ng(1:i-1)),ng(i),prod(ng(i+1:p)),na]; % bring arrays in vec shape
    shp = @(x) reshape(permute(reshape(x,sz),[2,1,3,4]),ng(i),[]);
    v = reshape(Mta,ng);
    for j=1:p, if i~=j, v = tmul(Kg.kron(j).factor,v,j); end, end
    if isnumeric(Kg.kron(i).factor)
      dhci = Kg.kron(i).dfactor( shp(v)*shp(Mtb)' );
    else
      kii = Kg.kron(i).factor.kii;
      [junk,ni] = apxGrid('expand',xg{i}); di = numel(ni);
      xs = cell(di,1);                % generic (Strang) circular embedding grid
      for j=1:di, n2 = floor(ni(j)-1/2)+1; xs{j} = [1:n2,n2-2*ni(j)+2:0]'; end
      [junk,dksi] = kii(apxGrid('expand',xs));
      Fvb = fftn(sum(fftn2(shp(v),ni',1).*fftn2(shp(Mtb), ni'),di+1));
      dhci = real(dksi(Fvb(:)));
    end
    dhyp = [dhyp; dhci(:)];
  end

function y = fftn2(x,n,t)             % fftn on trailing dimensions with padding
  if nargin<3, t = 0; end                                 % set a  default value
  nx = numel(x)/prod(n); y = reshape(x,[n(:)',nx]);                 % #instances
  if t
    for i=1:numel(n), y = fft(y,2*n(i)-1,i); end   % fftn on relevant dimensions
  else
    for i=1:numel(n), y = ifft(y,2*n(i)-1,i); end % ifftn on relevant dimensions
  end

% mode 11 (tmul)
% Multiplication with matrix/operator A along dimension dim of the tensor B.
% Note that tmul(A,B,1) = A*B if B is a matrix.
function C = tmul(A,B,dim)
if isnumeric(A), sa = size(A); else sa = A.size; end
sb = size(B); nb = ndims(B);
assert(dim>0 && floor(dim)==dim && dim<=nb)
assert(numel(sa)==2 && sb(dim)==sa(2))
if isnumeric(A), mvmA = @(x) A*x; else mvmA = A.mvm; end
if dim==1                              % along first dimension (save on permute)
  C = reshape(mvmA(reshape(B,sa(2),[])),[sa(1),sb(2:nb)]);
elseif dim==nb                          % along last dimension (save on permute)
  C = reshape(mvmA(reshape(B,[],sa(2))')',[sb(1:nb-1),sa(1)]);
else
  sb3 = [prod(sb(1:dim-1)),sa(2),prod(sb(dim+1:end))];
  C = permute(reshape(B,sb3),[2,3,1]);
  C = reshape(mvmA(reshape(C,sa(2),[])),[sa(1),sb3(3),sb3(1)]);
  C = reshape(permute(C,[3,1,2]),[sb(1:dim-1),sa(1),sb(dim+1:end)]);
end

% mode 12 (equi)
function eq = equi(xg,i)                        % grid along dim i is equispaced
  xi = xg{i};
  if iscell(xi)
    eq = true; for j=1:numel(xi), eq = eq && equi(xi,j); end
  else
    ni = size(xi,1);
    if ni>1                            % diagnose if data is linearly increasing
      dev = abs(diff(xi)-ones(ni-1,1)*(xi(2,:)-xi(1,:)));
      eq = max(dev(:))<1e-9;
    else
      eq = true;
    end
  end
  
% mode 13 (info)
function s = info(Kg,Mx,xg,deg)
  sk = 'K ='; Ks = Kg.kron; sk = [sk,' kron[ ']; p = numel(xg);
  for i=1:p
    if isnumeric(Ks(i).factor)
      si = sprintf('mat(%d)',size(Ks(i).factor,1));
    else
      sz = num2str(size(xg{i}{1},1));
      for j=2:numel(xg{i}), sz = [sz,'x',num2str(size(xg{i}{j},1))]; end
      si = sprintf('%s(%s)',Ks(i).factor.descr(1:4),sz);
    end
    if i<p, si = [si,' x ']; end
    sk = sprintf('%s%s',sk,si);
  end
  sk = sprintf('%s ]\n',sk);
  sm = 'M = ';
  xf = apxGrid('flatten',xg); deg = deg(:).*ones(numel(xf),1); id = 1;
  for i=1:p
    if apxGrid('equi',xg,i)
      si = sprintf('eq(d=%d',deg(id)); id = id+1;
      if iscell(xg{i})
        for j=2:numel(xg{i}), si=sprintf('%s,d=%d',si,deg(id)); id = id+1; end
      end
      si = [si,')'];
    else
      si = 'neq(d=1)'; id = id+1;
    end
    if i<p, si = [si,' x ']; end
    sm = sprintf('%s%s',sm,si);
  end
  sm = sprintf('%s, nnz=%d\n',sm,nnz(Mx));
  s = [sk,sm];

% mode 14 (Lanczos/ARPACK)
function [Q,T] = lanczos_arpack(B,v,d)     % perform Lanczos with at most d MVMs
  mv = ver('matlab');                                   % get the Matlab version
  if numel(mv)==0, error('No ARPACK reverse communication in Octave.'), end
  mv = str2double(mv.Version); old = mv<8;% Matlab 7.14=R2012a has old interface
  [junk,maxArraySize] = computer; m32 = maxArraySize==(2^31-1); % we have 32bit?
  if m32, intstr = 'int32'; else intstr = 'uint64'; end   % according data types
  intconvert = @(x) feval(intstr,x);  % init and allocate depending on # of bits
  n = length(v);
  v = bsxfun(@times,v,1./sqrt(sum(v.*v,1)));               % avoid call to normc
  ido = intconvert(0); nev = intconvert(ceil((d+1)/2)); ncv = intconvert(d+1);
  ldv = intconvert(n); info = intconvert(1);
  lworkl = intconvert(int32(ncv)*(int32(ncv)+8));
  iparam = zeros(11,1,intstr); ipntr = zeros(15,1,intstr);
  if exist('arpackc_reset'), arpackc_reset(); end
  iparam([1,3,7]) = [1,300,1]; tol = 1e-10;
  Q = zeros(n,ncv); workd = zeros(n,3); workl = zeros(lworkl,1); count = 0;
  
  while ido~=99 && count<=d
    count = count+1;
    if old
                   arpackc('dsaupd',ido,'I',intconvert(n),'LM',nev,tol,v,...
                                ncv,Q,ldv,iparam,ipntr,workd,workl,lworkl,info);
    else
      [ido,info] = arpackc('dsaupd',ido,'I',intconvert(n),'LM',nev,tol,v,...
                                ncv,Q,ldv,iparam,ipntr,workd,workl,lworkl,info);
    end
    if info<0
      error(message('ARPACKroutineError',aupdfun,full(double(info))));
    end
    if ido == 1, inds = double(ipntr(1:3));
    else         inds = double(ipntr(1:2)); end
    rows = mod(inds-1,n)+1; cols = (inds-rows)/n+1; % referenced column of ipntr
    if ~all(rows==1), error(message('ipntrMismatchWorkdColumn',n)); end
    switch ido                       % reverse communication interface of ARPACK
      case -1, workd(:,cols(2)) = B(workd(:,cols(1)));
      case  1, workd(:,cols(3)) =   workd(:,cols(1));
               workd(:,cols(2)) = B(workd(:,cols(1)));
      case  2, workd(:,cols(2)) =   workd(:,cols(1));
      case 99                                                        % converged
      otherwise, error(message('UnknownIdo'));
    end
  end
  ncv = int32(ncv);
  Q = Q(:,1:ncv-1);                                            % extract results
  T = diag(workl(ncv+1:2*ncv-1))+diag(workl(2:ncv-1),1)+diag(workl(2:ncv-1),-1);