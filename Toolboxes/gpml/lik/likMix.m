function [varargout] = likMix(lik, hyp, varargin)

% likMix - Mixture of likelihoods for regression/classification. 
% The expression for the likelihood is 
%   likMix(t) = sum_i=1..m  w_i * lik_i(t),
% where lik_i are the m individual likelihood functions combined by a weighted
% sum with weights wi.
%
% The hyperparameters are:
%
% hyp = [  log(w_1)
%          log(w_2)
%           ...
%          log(w_m)
%           hyp_1
%           hyp_2
%           ...
%           hyp_m ]
%
% Note that the mixture weights are normalised using the softmax function, so
% that 1 = sum w_i. Here, hyp_i are the hyperparameters for the individual
% likelihoods.
%
% Several modes are provided, for computing likelihoods, derivatives and moments
% respectively, see likFunctions.m for the details. In general, care is taken
% to avoid numerical issues when the arguments are extreme.
%
% Copyright (c) by Hannes Nickisch and Rowan McAllister, 2013-10-23.
%
% See also likFunctions.m.

m = numel(lik);                                   % number of mixture components
if m==0, error('We require at least one mixture component.'), end
for i = 1:m                                    % iterate over mixture components
  f = lik(i); if iscell(f{:}), f = f{:}; end    % expand cell array if necessary
  j(i) = cellstr(feval(f{:}));                           % collect number hypers
end

s = @(x) logsumexp2(x);                 % define a shortcut to util/logsumexp2.m

if nargin<4                                        % report number of parameters
  L = char(j(1)); for i=2:length(lik), L = [L, '+', char(j(i))]; end
  varargout = {['(',num2str(m),'+',L,')']}; return
end
mu = varargin{2}; n = length(mu);                       % number of query points
if numel(varargin)>3
  inf = varargin{4};                            % varargin = {y, mu, s2, inf, i}
else
  inf = '';
end
v = [];               % v vector indicates to which likelihood parameters belong
nhyp = m;       % number of hyperparameters; the first m are the mixture weights
for i = 1:m
  ni = eval(char(j(i))); nhyp = nhyp + ni; v = [v, repmat(i,1,ni)];
end
if nhyp>length(hyp), error('not enough hyperparameters'), end
if strcmp(inf,'infVB'), error('infVB not supported'); end

na = nargout;
if nargout>1 && strcmp(inf,'infLaplace') && nargin>6, na = max(na,3); end
argout = zeros(n,m,na);  % collect output arguments from individual lik fun
for i=1:m                                    % iteration over mixture components
  f = lik(i); if iscell(f{:}), f = f{:}; end    % expand cell array if necessary
  id = find(v==i)+m;                                    % hyperparameter indices
  out = cell(size(argout,3),1); nn = min(nargin,6)-2;    % do not ask for derivs
  [out{:}] = feval(f{:}, hyp(id), varargin{1:nn});  % ask individual likelihoods
  for j=1:size(argout,3), argout(:,i,j) = out{j}; end
end

lw = reshape(hyp(1:m),1,[]); lwn = lw-s(lw);   % w = exp(lwn), sum(w)=1, weights
one_lwn = ones(n,1)*lwn;  % the same as above only that the weights are repeated
varargout = cell(nargout,1);              % allocate memory for output arguments
if nargin<6                              % prediction mode if inf is not present
  % lp and ymu outputs are of the weighted sum of integral structure
  if nargout>0, varargout{1} = s(one_lwn + argout(:,:,1)); end
  if nargout>1, varargout{2} = sum(exp(one_lwn).*argout(:,:,2), 2); end
  if nargout>2
    % ys2 output computed using the parallel axis theorem
    argout3corr = argout(:,:,3) + (argout(:,:,2)-varargout{2}*ones(1,m)).^2;
    varargout{3} = sum(exp(one_lwn).*argout3corr, 2);
  end
else
  if nargin>6, i = varargin{5}; end                   % make arguments available
  switch inf 
  case {'infEP','infLaplace'}           % the two cases are structurally similar
    if nargin<7                                             % no derivative mode
      if nargout>0                   % lf={lZ,lp} for inf={'infEP','infLaplace'}
        lfi = argout(:,:,1);     % partition functions of individual likelihoods
        lf = s(lfi + one_lwn);                  % weighted sum in the exp domain
        varargout{1} = lf;
        if nargout>1                                                       % dlf
          dlfi = argout(:,:,2);                % derivatives of individual terms
          % Using f*dlf=df, f=exp(lf), f=sum_i wi*fi, df=sum_i wi*dfi we obtain
          %   dlf = sum_i exp(lfi-lf+lwi)*dlfi = sum_i ai*dlfi.
          a = exp(lfi - lf*ones(1,m) + one_lwn);
          dlf = sum(a.*dlfi,2);
          varargout{2} = dlf;
          if nargout>2                                                    % d2lf
            d2lfi = argout(:,:,3);         % 2nd derivatives of individual terms
            % Using d2lf=d(df/f)=(d2f*f-df^2)/f^2 <=> d2f=f*(d2lf+dlf^2) and
            % d2f = sum_i wi*d2fi, we get d2lf = sum_i ai*(d2lfi+dlfi^2)-dlf^2.
            d2lf = sum(a.*(d2lfi+dlfi.*dlfi),2) - dlf.*dlf;
            varargout{3} = d2lf;
            if nargout>3 && strcmp(inf,'infLaplace')                      % d3lf
              % Using d3f = sum_i wi*d3fi, d3lf=d((d2f*f-df^2)/f^2), we obtain
              % d3f = f*(d3lf+dlf*(3*d2lf+dlf^2)) <=>
              % d3lf = d3f/f - dlf*(3*d2lf+dlf^2).
              d3lfi = argout(:,:,4);       % 3rd derivatives of individual terms
              d3lf = sum(a.*(d3lfi+dlfi.*(3*d2lfi+dlfi.*dlfi)),2);
              d3lf = d3lf - dlf.*(3*d2lf+dlf.*dlf);
              varargout{4} = d3lf;
            end
          end
        end
      end
    else                                           % d(lf)/dhyp and d(d2lf)/dhyp
      if nargout>0
        lfi = argout(:,:,1);     % partition functions of individual likelihoods
        lf = s(lfi + one_lwn);                  % weighted sum in the exp domain
        a = exp(lfi - lf*ones(1,m) + one_lwn);
        if i<=m                % the first m hyperparameters are mixture weights
          lf_dhyp = a(:,i) - exp(lwn(i));
          varargout{1} = lf_dhyp;
          if nargout>1 && strcmp(inf,'infLaplace')                 % d(dlf)/dhyp
            dlfi = argout(:,:,2);
            dlwn = zeros(size(lw)); dlwn(i) = 1;
            dlwn = dlwn - exp(lw(i)-s(lw));
            dla = ones(n,1)*dlwn - lf_dhyp*ones(1,m);          % d log(v) / dhyp
            da = a.*dla;
            dlf_dhyp = sum(da.*dlfi,2);
            varargout{2} = dlf_dhyp;
            if nargout>2                                          % d(d2lf)/dhyp
              d2lfi = argout(:,:,3);
              dlf = sum(a.*dlfi,2);
              d2lf_dhyp = sum(da.*(d2lfi+dlfi.*dlfi),2) - 2*dlf.*dlf_dhyp;
              varargout{3} = d2lf_dhyp;
            end
          end
        else
          j = v(i-m);           % index of mixture component under consideration
          id = find(v==j)+m;                           % hyper parameter indices
          i = i-min(id)+1;          % index of likelihood parameter of mixture j
          f = lik(j); if iscell(f{:}), f = f{:}; end  % expand cell if necessary
          if nargout==1
            lfi_dhyp = feval(f{:}, hyp(id), varargin{1:end-1}, i);
          elseif nargout==2
            [lfi_dhyp,dlfi_dhyp] = feval(f{:}, hyp(id), varargin{1:end-1}, i);
          else
            [lfi_dhyp,dlfi_dhyp,d2lfi_dhyp] = ...
                                     feval(f{:}, hyp(id), varargin{1:end-1}, i);
          end
          lf_dhyp = exp(lwn(j)+lfi(:,j)-lf) .* lfi_dhyp(:);         % d(lf)/dhyp
          varargout{1} = lf_dhyp;
          if nargout>1                                             % d(dlf)/dhyp
            dlfi = argout(:,:,2); d2lfi = argout(:,:,3);
            dlf = sum(a.*dlfi,2);
            dla = -lf_dhyp*ones(1,m); dla(:,j) = dla(:,j)+lfi_dhyp(:); 
            da = a.*dla; 
            dlf_dhyp = sum(da.*dlfi,2) + a(:,j).*dlfi_dhyp(:);
            varargout{2} = dlf_dhyp;
            if nargout>2                                          % d(d2lf)/dhyp
              p = d2lfi+dlfi.*dlfi;
              dpj = 2*dlfi(:,j).*dlfi_dhyp(:) + d2lfi_dhyp(:);
              d2lf_dhyp = a(:,j).*dpj + sum(da.*p,2) - 2*dlf.*dlf_dhyp;
              varargout{3} = d2lf_dhyp;
            end
          end
        end
      end
    end

  case 'infVB'
    error('infVB not supported')
  end
end