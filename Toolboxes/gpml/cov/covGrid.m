function [K,Mx,xe] = covGrid(cov, xg, hyp, x, z, i)

% covGrid - Kronecker covariance function based on a grid.
%
% The grid g is represented by its p axes xg = {x1,x2,..xp}. An axis xi is of
% size (ni,di) and the grid g has size (n1,n2,..,np,D), where D=d1+d2+..+dp.
% Hence, the grid contains N=n1*n2*..*np data points. The axes do neither need
% to be sorted in any way nor do they need to be 1d i.e. ni>=1.
%
% The covGrid function can be used to expand the cell array xg into an expanded
% multivariate grid xe of size (N,D) via:
%     [xe,nx,Dx]  = covGrid('expand',xg);
% The operation can be reverted by:
%     xg = covGrid('factor',{xe,ng,Dg});
%
% The variables v={x,z} can either be a) grid indices or b) data points.
% a) The variable v has size (nv,1) and contains integers from [1,N]. Then
%    the datapoints are obtained as g2 = reshape(g,N,D); v = g2(v,:).
% b) The variable v has size (nv,D) and directly represents the data points.
% The mechanism works for x and z separately.
% 
% The resulting covariance matrix is given by:
%   K = kron( kron(...,K{2}), K{1} ) = K_p x .. x K_2 x K_1.
%
% The hyperparameters are:
% hyp = [ hyp_1
%         hyp_2
%          ..
%         hyp_p ],
%
% Copyright (c) by Hannes Nickisch and Andrew Wilson 2014-12-04.
%
% See also COVFUNCTIONS.M, INFGRID.M.

if nargin<2, error('Not enough parameters provided.'), end
if     strcmp(cov,'expand')  % convert between full grid and axes representation
  [K,Mx,xe] = expandgrid(xg); return
elseif strcmp(cov,'factor')
  K = factorgrid(xg{:}); return
end

p = numel(xg); ng = zeros(p,1); Dg = zeros(p,1);   % number of Kronecker factors
if numel(cov)~=p, error('We require p factors.'), end                               
for ii = 1:p                                 % iterate over covariance functions
  [ng(ii),Dg(ii)] = size(xg{ii});
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
if nargin==6 && i>length(v), error('Unknown hyperparameter'), end

toep = true(p,1);     % grid along dimension is equispaced -> Toeplitz structure
for ii=1:p
  if ng(ii)>1       % diagnose Toeplitz structure if data is linearly increasing
    dev = abs(diff(xg{ii})-ones(ng(ii)-1,1)*(xg{ii}(2,:)-xg{ii}(1,:)));
    toep(ii) = max(dev(:))<1e-9;
  end
end

N = prod(ng); n = size(x,1); D = sum(Dg);     % expanded grid and data dimension
ix = isidx(x,N);               % determine whether x is an index or a data array
if dg               % evaluate as full dense vector for diagonal covariance case
  K = 1;                       % xg is not assumed to form a grid for z = 'diag'
  for ii = 1:length(cov)                       % iteration over factor functions
    f = cov(ii); if iscell(f{:}), f = f{:}; end % expand cell array if necessary
    d = sum(Dg(1:ii-1))+(1:Dg(ii));                     % dimensions of interest
    if nargin<6, i = 0; vi = 0; else vi = v(i); end; % which covariance function
    if i<=length(v)
      if ix, xii = xg{ii}; else xii = x(:,d); end  % switch Kronecker/plain prod
      if ii==vi
        j = sum(v(1:i)==vi);                % which parameter in that covariance
        Kj = feval(f{:}, hyp(v==ii), xii, z, j);        % deriv Kronecker factor
      else
        Kj = feval(f{:}, hyp(v==ii), xii, z);           % plain Kronecker factor
      end
      if ix, K = kron(K,Kj); else K = K.*Kj; end   % switch Kronecker/plain prod
    else error('Unknown hyperparameter')
    end
  end
  if ix, K = K(x); end, return
end
if ~ix, error('Off-grid input only supported for z.'), end

if isidx(z,N); Mz = z; z = covGrid('expand',xg); z = z(Mz,:); end
K = cell(p,1);                                    % covariance Kronecker factors
for ii = 1:length(cov)                         % iteration over factor functions
  f = cov(ii); if iscell(f{:}), f = f{:}; end   % expand cell array if necessary
  d = sum(Dg(1:ii-1))+(1:Dg(ii));                       % dimensions of interest
  if isnumeric(z) && ~isempty(z)                                   % cross terms
    zd = z(:,d);
  elseif xeqz && toep(ii)                           % we have Toeplitz structure
    zd = xg{ii}(1,:);
  else                                                        % symmetric matrix
    zd = z;
  end
  if nargin<6, i = 0; vi = 0; else vi = v(i); end;   % which covariance function
  if i<=length(v)
    if ii==vi
      j = sum(v(1:i)==vi);                  % which parameter in that covariance
      K{ii} = feval(f{:}, hyp(v==ii), xg{ii}, zd, j);   % deriv Kronecker factor
    else
      K{ii} = feval(f{:}, hyp(v==ii), xg{ii}, zd);      % plain Kronecker factor
    end
  else error('Unknown hyperparameter')
  end
  if xeqz && toep(ii), K{ii} = {'toep',K{ii}}; end      % make Toeplitz explicit
end

if ~xeqz                                                    % expand cross terms
  Ks = K; K = Ks{1}; for ii = 2:p, K = kron1(Ks{ii},K); end
  if ix && (numel(x)~=N || max(abs(x-(1:N)'))>0), K = K(x,:); end
end
if nargout>1, if ix, Mx = sparse(1:n,x,1,n,N); end, end
if nargout>2, xe = covGrid('expand',xg); end

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

function [x,ng,Dg] = expandgrid(xg)                    % expand a Kronecker grid
  p = numel(xg); x = xg{1};                                 % expanded grid data
  ng = zeros(p,1); Dg = zeros(p,1); [ng(1),Dg(1)] = size(xg{1});
  for i=2:p
    szx = size(x); [ng(i),Dg(i)] = size(xg{i});
    xold = repmat(reshape(x,[],1,szx(end)),[1,ng(i),1]);
    xnew = repmat(reshape(xg{i},[1,ng(i),Dg(i)]),[size(xold,1),1,1]);
    x = reshape(cat(3,xold,xnew),[szx(1:end-1),ng(i),szx(end)+Dg(i)]);
  end
  x = reshape(x,[],size(x,ndims(x)));

function xg = factorgrid(x,ng,Dg)                      % factor a Kronecker grid
 p = numel(ng); xg = cell(p,1);          % extract individual grid components xg
for i=1:p
  x = reshape(x,[prod(ng(1:i-1)), ng(i), prod(ng(i+1:end)), sum(Dg)]);
  xg{i} = reshape(x(1,:,1,sum(Dg(1:i-1))+(1:Dg(i))), ng(i), Dg(i));
end