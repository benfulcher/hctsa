function out = sgnchange(y)
% Returns indices where the input vector, y, changes sign
% Ben Fulcher
    
out = find(y(2:end).*y(1:end-1) < 0);

end