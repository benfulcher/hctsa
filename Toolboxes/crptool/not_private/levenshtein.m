function y=levenshtein(a,b)
% LEVENSHTEIN   Calculates the Levenshtein distance.
%     D=LEVENSHTEIN(A,B) calculates the Levenshtein distance
%     between row vectors A and B. If A and B are arrays,
%     Levenshtein distances are calculated for all rows in
%     these array, thus D is a column vector of length(A).
%     The size(A,1) has to match size(B,1).
%
%     Reference:
%     Levenshtein, V. I.: 
%     Binary codes capable of correcting deletions, insertions, 
%     and reversals, Doklady Akademii Nauk SSSR, 163(4), 1965, 
%     845-848, 1965 (Russian). 
%     English in: Soviet Physics Doklady, 10(8), 1966, 707-710.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:36:04 $
% $Revision: 4.2 $

if size(a,1) ~= size(b,1)
    error('Array''s dimension 1 does not match.')
end


Na = size(a,2); Nb = size(b,2);

y = zeros(size(a,1),1);

for k = 1:size(a,1)
    D = zeros([Na+1, Nb+1]);
    D(1,:) = 0:Nb; D(:,1) = 0:Na;

    for i = 2:Na + 1
       for j = 2:Nb + 1
          cost = (b(k, j-1) ~= a(k, i-1));
          D(i,j)=min([D(i-1,j-1) + cost, D(i-1,j)+1, D(i,j-1)+1]);
       end
    end
    y(k) = D(end,end);
end
