function r = KP_inform_ch(y,nbins)
% INFORM calculates the information given a probability distribution
% in the form of bincounts.  
% Example usage:  (s,t) = hist(data); inform(s)
% r is normalized to the information in the uniform distribution
% with length(bincnts)
% Copyright (c) 1996 by D. Kaplan, All Rights Reserved

[bincnts dnx] = hist(y,nbins);

cnts = bincnts(bincnts>0);
% total = sum(cnts);
cnts = cnts/sum(cnts);
r = - sum(cnts.*log(cnts))/log(length(bincnts));