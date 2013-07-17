function y = minimize_lbfgsb_objfun(X,varargin)

  % extract input arguments
  varargin = varargin{1}; strctX = varargin{2}; f = varargin{1};

  % global variables serve as communication interface between calls
  global minimize_lbfgsb_iteration_number
  global minimize_lbfgsb_objective
  global minimize_lbfgsb_gradient
  global minimize_lbfgsb_X

  recompute = 0;
  if minimize_lbfgsb_iteration_number==0
    recompute = 1;
  else
    if norm(unwrap(X)-unwrap(minimize_lbfgsb_X))>1e-10
      recompute = 1;
    end
  end
  if recompute
    [y,G] = feval(f,rewrap(strctX,X),varargin{3:end});
  else
    y = minimize_lbfgsb_objective(minimize_lbfgsb_iteration_number);
    G = minimize_lbfgsb_gradient;
  end
  
  % increase global counter, memorise objective function, gradient and position
  minimize_lbfgsb_iteration_number = minimize_lbfgsb_iteration_number+1;
  minimize_lbfgsb_objective(minimize_lbfgsb_iteration_number,1) = y;
  minimize_lbfgsb_gradient = G;
  minimize_lbfgsb_X = X;


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
