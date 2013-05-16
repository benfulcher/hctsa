function cout = intermutual(cin1, cin2, n)

%tstoolbox/@core/intermutual
%   Syntax:
%     * intermutual(cin1,cin2,n)
%
%   Input Arguments:
%     * cin1,cin2 - core objects
%
%   Calculates the mutual information of cin1 and cin2.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


N = dlens(cin1,1);
cout = zeros(n,1);
epsilon = 1e-9;
points1 = data(cin1);
points1 = points1 - min(points1);
points2 = data(cin2);
points2 = points2 - min(points2);
ma1 = max(points1);
ma2 = max(points2);
increment = 1/N;

for i=1:n
	partitionen = 2^i;
	pointsA = 1+floor(points1 / (ma1/(partitionen-epsilon)));
	pointsB = 1+floor(points2 / (ma2/(partitionen-epsilon)));
	% Aus den Werten sind jetzt fertige Indizes fuer die Partititonen geworden
	
	histA = sparse(pointsA, ones(N,1), increment);
	histB = sparse(ones(N,1), pointsB, increment);
	histAB = sparse(pointsA, pointsB, increment);
	[ind1, ind2, value] = find(histAB);
	amf = sum(value .* log2(value ./ (histA(ind1) .* histB(ind2)')));	% '

	cout(i) = amf;
end

%cout = core(amf);

	
