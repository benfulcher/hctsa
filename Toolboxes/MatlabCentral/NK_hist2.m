function [MN, xedges, yedges] = NK_hist2(x, y, xedges, yedges)
% function MN = NK_hist2(x, y, xedges, yedges)
%
% 2D histogram: Extract the number of joint events - (x,y) data value pairs 
% that fall in each bin of the grid defined by xedges and yedges.
%
% ====================================================================
% Uses:
%    [N,BIN] = histc(X,EDGES)
% which returns
%   1) N is a LENGTH(EDGES) vector, N(k) will count the value X(i) if EDGES(k) <= X(i) < EDGES(k+1).  The
%    last bin will count any values of X that match EDGES(end).
%   2) an index matrix BIN.
%    If X is a vector, N(K) = SUM(BIN==K).
%    BIN is zero for out of range values.
%    If X is an m-by-n matrix, then, for j=1:n, N(K,j) = SUM(BIN(:,j)==K); end
% ====================================================================
%
%      Please, see the notes to histc too.
%      N.B. It is always a better idea to use the
%      HISTC mex (a much faster compiled C code) if you have it
%      Then just replace the histc with HISTC in all calls
%      contained in the NK_hist2() .m function
%
% (c) Nedialko Krouchev 2006, Universite de Montreal, GRSNC

if nargin ~= 4
    error('The four input arguments are required!');
end
if any(size(x) ~= size(y))
    error('The size of the two first input vectors should be same!');
end

[xn, xbin] = histc(x,xedges);
[yn, ybin] = histc(y,yedges);

xnbin = length(xedges);
ynbin = length(yedges);

% xbin, ybin are zero for out of range values
kkL = find( x<xedges(1) ); kkR = find( x>xedges(xnbin) );
if ~isempty( kkL ),
xbin = xbin + 1; xnbin = xnbin + 1;
xedges(2:xnbin) = xedges; xedges(1) = min(x);
end
xbin( kkL ) = 1;
if ~isempty( kkR ),
xnbin = xnbin + 1; xedges(xnbin) = max(x);
end
xbin( kkR ) = xnbin;

kkL = find( y<yedges(1) ); kkR = find( y>yedges(ynbin) );
if ~isempty( kkL ), ybin = ybin + 1; ynbin = ynbin + 1;
yedges(2:ynbin) = yedges; yedges(1) = min(y);
end
ybin( kkL ) = 1;
if ~isempty( kkR ), ynbin = ynbin + 1; yedges(ynbin) = max(y); end
ybin( kkR ) = ynbin;

xyBinEdges = 1:xnbin*ynbin;

% ====================================================================
% A more Elegant end-spiel:
%
% If x  belongs to jBin=xbin(x), and y  belongs to iBin=ybin(y),
% Then (x,y) pairs belong "columnwise" to:
%       ijBin = (jBin-1)*ynbin + iBin

xyBin = (xbin-1)*ynbin + ybin;


% ====================================================================

MN = histc(xyBin, xyBinEdges);
MN = reshape(MN, ynbin, xnbin);
end