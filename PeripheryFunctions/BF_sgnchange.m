function out = BF_sgnchange(y,dofind)
% Returns indices where the input vector, y, changes sign
% Ben Fulcher

if nargin < 2
    dofind = 0; % by default don't find, just return logical of where the sign changes (faster)
end

out = (y(2:end).*y(1:end-1) < 0); % successive differences change sign

% [case of adding equality, as <= 0, but I think better to 
% make hard changes in sign as < 0]

if dofind
    out = find(out); % return indicies of where sign changes
end

end