function out=ST_mm(y,s)
    switch s
        case 'max'
            out=max(y);
        case 'min'
            out=min(y);
    end
end