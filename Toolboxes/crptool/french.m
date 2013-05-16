function h = french(m,flag)
%FRENCH    French's flag color map.
%   FRENCH(M) returns an M-by-3 matrix containing a "french flag" colormap.
%   FRENCH, by itself, is the same length as the current colormap.
%
%   FRENCH(M,FLAG) enables a modification of the colormap, where
%   FLAG can be either 1 (default) or 2.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(french)
%
%   See also HOT, HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG, 
%   COLORMAP, RGBPLOT.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2002-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:32:57 $
% $Revision: 2.3 $
%
% $Log: french.m,v $
% Revision 2.3  2009/03/24 08:32:57  marwan
% copyright address changed
%
% Revision 2.2  2006/03/16 13:52:43  marwan
% code refreshed
%
% Revision 2.1  2004/11/10 07:07:51  marwan
% initial import
%


if nargin < 1
    m = size(get(gcf,'colormap'),1); 
else
    if isempty(m)
        m = size(get(gcf,'colormap'),1); 
    end
end
if nargin < 2, flag = 1; end

n1 = fix(3*m/8);
n2 = fix(m/4);
n3 = fix(m/2);

switch flag
    case 1
        r = [ones(n3,1); (sqrt(1-((1:n3)/n3).^2))'];
        g = [flipud( (sqrt(1-((1:n3)/n3).^2))'); (sqrt(1-((1:n3)/n3).^2))'];
        b = [flipud((sqrt(1-((1:n3)/n3).^2))'); ones(n3,1)];
    case 2
        r = [ones(n1+n2,1); (n1-1:-1:0)'/n1;];
        g = [(0:n1-1)'/n1; ones(n2,1); (n1-1:-1:0)'/n1;];
        b = [(0:n1-1)'/n1; ones(n1+n2,1);];
end

h = [r g b];

if size(h,1) < m
    h(ceil(m/2)+1:m,:) = h(ceil(m/2):end,:);
    h(ceil(m/2),:) = 1;
end
