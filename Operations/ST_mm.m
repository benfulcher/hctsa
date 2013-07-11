function out = ST_mm(y,s)
% Ben Fulcher, 2008

switch s
    case 'max'
        out = max(y);
    case 'min'
        out = min(y);
    otherwise
        error('Unknown method ''%s''',s)
end

end