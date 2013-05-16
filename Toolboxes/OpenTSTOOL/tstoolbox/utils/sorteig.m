function [v,svals] = sorteig(C, percent)

% same as eig, except eigenvectors in sorted order, highest eigenvalue first
% and eigenvalue are just a vector, not a matrix
% and only some eigenvectors (depending on percent) are computed

[v,d] = eig(C);
vals = diag(d);
[svals, index] = sort(vals(:)); 
svals = flipud(svals);
index = flipud(index);


rlvm = min(find(cumsum(svals) >= (sum(svals)*percent/100))); %percent; 

if isempty(rlvm)			% in case percentage was choosen over 100 %,
	rlvm = length(svals); % return all eigenvalues
end


svals = svals(1:rlvm);
v = v(:, index(1:rlvm));

%svals
%bar(100*svals/sum(svals))
%v
