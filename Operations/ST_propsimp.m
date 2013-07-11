function out = ST_propsimp(x,anuey)

N=length(x);

switch anuey
    case 'zeros' % returns the proportion of zeros in the input vector
        out=length(find(x == 0))/N;
    case 'positive'
        out=length(find(x>0))/N;
    case 'geq0'
        out=length(find(x>=0))/N;
end


end