function [X, fX, i] = minimize_lbfgsb(X, f, length, varargin)

% Minimize a differentiable multivariate function using quasi Newton.
%
% Usage: [X, fX, i] = minimize_lbfgsb(X, f, length, P1, P2, P3, ... )
% 
% X       initial guess; may be of any type, including struct and cell array
% f       the name or pointer to the function to be minimized. The function
%         f must return two arguments, the value of the function, and it's
%         partial derivatives wrt the elements of X. The partial derivative  
%         must have the same type as X.
% length  length of the run; if it is positive, it gives the maximum number of
%         line searches, if negative its absolute gives the maximum allowed
%         number of function evaluations. Optionally, length can have a second
%         component, which will indicate the reduction in function value to be
%         expected in the first line-search (defaults to 1.0).
% P1, P2  ... parameters are passed to the function f.
%
% X       the returned solution
% fX      vector of function values indicating progress made
% i       number of iterations (line searches or function evaluations, 
%         depending on the sign of "length") used at termination.
%
% The function returns when either its length is up, or if no further progress
% can be made (ie, we are at a (local) minimum, or so close that due to
% numerical problems, we cannot get any closer). NOTE: If the function
% terminates within a few iterations, it could be an indication that the
% function values and derivatives are not consistent (ie, there may be a bug in
% the implementation of your "f" function).
%
% Copyright (C) 2010 by Hannes Nickisch, 2010-02-05

% global variables serve as communication interface between calls
global minimize_lbfgsb_iteration_number
global minimize_lbfgsb_objective
global minimize_lbfgsb_gradient
global minimize_lbfgsb_X

% init global variables
minimize_lbfgsb_iteration_number = 0;
minimize_lbfgsb_objective = Inf;
minimize_lbfgsb_gradient = 0*unwrap(X);
minimize_lbfgsb_X = X;

X0 = X;
lb = -Inf*ones(size(unwrap(X0)));
ub =  Inf*ones(size(unwrap(X0)));
maxiter = abs(length); % max number of iterations

% no callback routine used so far
% m is the number of saved vectors used to estimate the Hessian
% factr is the precision 1e-12
X = lbfgsb( unwrap(X0), lb, ub, 'minimize_lbfgsb_objfun', ...
                                'minimize_lbfgsb_gradfun', ...
         {f,X0,varargin{:}}, [], ...
         'maxiter',maxiter, 'm',4, 'factr',1e-12, 'pgtol',1e-5);
i = minimize_lbfgsb_iteration_number;
fX = minimize_lbfgsb_objective;
X = rewrap(X0,X);

% clear global variables
clear minimize_lbfgsb_iteration_number
clear minimize_lbfgsb_objective
clear minimize_lbfgsb_gradient
clear minimize_lbfgsb_X


% Extract the numerical values from "s" into the column vector "v". The
% variable "s" can be of any type, including struct and cell array.
% Non-numerical elements are ignored. See also the reverse rewrap.m. 
function v = unwrap(s)
  v = [];   
  if isnumeric(s)
    v = s(:);                       % numeric values are recast to column vector
  elseif isstruct(s)
    v = unwrap(struct2cell(orderfields(s)));% alphabetize, conv to cell, recurse
  elseif iscell(s)
    for i = 1:numel(s)            % cell array elements are handled sequentially
      v = [v; unwrap(s{i})];
    end
  end                                                  % other types are ignored

% Map the numerical elements in the vector "v" onto the variables "s" which can
% be of any type. The number of numerical elements must match; on exit "v"
% should be empty. Non-numerical entries are just copied. See also unwrap.m.
function [s v] = rewrap(s, v)
  if isnumeric(s)
    if numel(v) < numel(s)
      error('The vector for conversion contains too few elements')
    end
    s = reshape(v(1:numel(s)), size(s));           % numeric values are reshaped
    v = v(numel(s)+1:end);                       % remaining arguments passed on
  elseif isstruct(s) 
    [s p] = orderfields(s); p(p) = 1:numel(p);     % alphabetize, store ordering  
    [t v] = rewrap(struct2cell(s), v);                % convert to cell, recurse
    s = orderfields(cell2struct(t,fieldnames(s),1),p); % conv to struct, reorder
  elseif iscell(s)
    for i = 1:numel(s)            % cell array elements are handled sequentially 
      [s{i} v] = rewrap(s{i}, v);
    end
  end                                            % other types are not processed
