function r = inform(bincnts)
% INFORM calculates the information given a probability distribution
% in the form of bincounts.  
% Example usage:  (s,t) = hist(data); inform(s)
% r is normalized to the information in the uniform distribution
% with length(bincnts)
% Copyright (c) 1996 by D. Kaplan, All Rights Reserved

g = bincnts > 0;
cnts = bincnts(find(g));
total = sum(cnts);
cnts = cnts/total;
r = - sum( cnts.*log(cnts))/log(length(bincnts));



