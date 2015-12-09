function [f,g] = NetVis_Jtilda2d(x,y,n)
% does estimation of the repulsion forces
% using RANK and ROTATE algorithm
% Code provided by David Smith

N = length(x);
if(nargin<3)
    n=10;
end
multfactor = pi/(n*4);

% Initialization:
xrank = zeros(N,1);
yrank = zeros(N,1);

p = pi/(2*n);
[~, xind] = sort(x);
xrank(xind) = 1:N;
f = 2*xrank-(N+1);
[~, yind] = sort(y);
yrank(yind) = 1:N;
g = 2*yrank-(N+1);

for i = 1:n-1
    c = cos(i*p);
    s = sin(i*p);
    x2 = c*x-s*y;
    y2 = s*x+c*y;

    [~, xind] = sort(x2);
    xrank(xind) = 1:N;
    f1 = 2*xrank-(N+1);
    [~, yind] = sort(y2);
    yrank(yind) = 1:N;
    g1 = 2*yrank-(N+1);

    f2 = c*f1+s*g1;
    g2 = -s*f1+c*g1;
    f = f+f2;
    g = g+g2;
end

f = f*multfactor;
g = g*multfactor;

end
