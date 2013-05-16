function y=dtw(a,b)
% DTW   Calculates the dynamic time warping distance.
%     D=DTW(A,B) calculates the dynamic time warping (DTW) 
%     distance between row vectors A and B. If A and B are 
%     arrays, DTW distances are calculated for all rows in
%     these array, thus D is a column vector of length(A).
%     The size(A,1) has to match size(B,1).
%
%     Reference:
%     Myers, C. S, Rabiner, L. R.:
%     A comparative study of several dynamic time-warping algorithms 
%     for connected word recognition, The Bell System Technical 
%     Journal, 60(7), 1982, 1389-1409.
%     FastDTW: Toward Accurate Dynamic Time Warping in Linear Time 
%     and Space, Stan Salvador and Philip Chan. Intelligent 
%     Data Analysis, 2007.

% Copyright (c) 2008 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2008/02/27 11:21:22 $
% $Revision: 4.1 $


if size(a,1) ~= size(b,1)
    error('Array''s dimension 1 does not match.')
end


Na = size(a,2); Nb = size(b,2);

y = zeros(size(a,1),1);
%h = waitbar(0, 'Calculation DTW distance');
for k = 1:size(a,1), %if k/100 == fix(k/100), waitbar(k/size(a,1)), end
    aa = a(k,:)'; bb = b(k,:)';
    eucDis = (repmat(aa,1,Nb) - repmat(bb',Na,1)).^2; 

    D = zeros(size(eucDis));
    D(1,1) = eucDis(1,1);
    D(2:Na,1) = eucDis(2:Na,1) + D((2:Na)-1,1);
    D(1,2:Nb) = eucDis(1,2:Nb) + D(1,(2:Nb)-1);
    for n = 2:Na
        for m = 2:Nb
            D(n,m) = eucDis(n,m) + min([D(n-1,m), D(n-1,m-1), D(n,m-1)]);
        end
    end

    y(k) = D(Na,Nb);
end
%delete(h)
