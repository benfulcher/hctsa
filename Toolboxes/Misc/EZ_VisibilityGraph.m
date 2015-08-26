function [A] = visibilitygraph(T, Nmax)
% input:
% T: time series
% Nmax: the size of the time series T
% output:
% A: Adjacent matrix of the corresponding visibility graph
%
% Enyu ZHUANG (Zoey) 23/9/13

A = zeros(Nmax, Nmax);

for i = 1:(Nmax-1)
    for j = i+1:Nmax
        if(j == i+1)
            A(i,j)=1;
            A(j,i)=1;
            S(i, j) = T(j)-T(i);
        else
            S(i, j) = (T(j)-T(i))/(j-i);
            for k = i+1: j-1
                if(S(i,k) >= S(i,j))
                    break;
                elseif(k == (j-1))
                    A(i,j)=1;
                    A(j,i)=1;
                end
            end
        end
    end
end
