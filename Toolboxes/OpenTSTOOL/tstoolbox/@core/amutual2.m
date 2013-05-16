function cout = amutual(cin, len)

%tstoolbox/@core/amutual2
%   Syntax:
%     * amutual(cin, len)
%
%   Input Arguments:
%     * cin core object
%     * len maximal lag
%
%   amutual2 calculates the mutual information of a time series against
%   itself, with increasing lag uses equidistant partitioning to compute
%   histograms.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

partitionen = 128;
epsilon = 1e-10;

N = dlens(cin,1);
points = data(cin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cmerk 23.3.1998
% Transformation auf Rangzahlenfolge



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
points = points - min(points);
points = 1+floor(points / (max(points)/(partitionen-epsilon)));
% Aus den Werten sind jetzt fertige Indizes fuer die Partititonen geworden

amf = zeros(len,1);
overlap = N+1-len;
increment = 1/overlap;
one = ones(overlap,1);
histA= sparse(points(1:overlap), one, increment);

for lag=0:len-1	
	histB = sparse(one, points(1+lag:overlap+lag), increment);
	histAB = sparse(points(1:overlap), points(1+lag:overlap+lag), increment);
	[ind1, ind2, value] = find(histAB);
	amf(lag+1) = sum(value .* log2(value ./ (histA(ind1) .* histB(ind2)')));	% '
end

cout = core(amf);

	
% 	s = 0;		
% 	for i=1:length(ind1)
% 		s = s + value(i) * log2(value(i)/(histA(ind1(i))*histB(ind2(i))));
% 	end
	
% 	sum(histA)  
% 	sum(histB)
% 	sum(sum(histAB))

% 	sp2 = sparse(histA * histB);
% 	index = find(sp2);
% 	quot = histAB;
% 	quot(index) = quot(index) ./  sp2(index);	% in dieser Matrix stehen die Quotienten aus P(AB)/(P(A)*P(B))
% 									
% 	index = find(quot); % Nur Elemente logarithmieren, die ungleich Null sind
% 	quot(index) = log2(quot(index));	
% 
% 	amf(lag+1) = sum(sum(histAB(index) .* quot(index)));	% Gewichtung und Aufsummierung
% 	amf(lag+1) = s;
