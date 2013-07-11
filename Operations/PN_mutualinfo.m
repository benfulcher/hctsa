function h = PN_mutualinfo(vec1,vec2)
% == == == == == == == == == == == == == == == == == == == == == == == == == == == ===
%
%This is a prog in the MutualInfo 0.9 package written by 
% Hanchuan Peng.
%
%Disclaimer: The author of program is Hanchuan Peng
%      at <penghanchuan@yahoo.com> and <phc@cbmv.jhu.edu>.
%
%The CopyRight is reserved by the author.
%
%Last modification: April/19/2002
%
% == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
%
% h = mutualinfo(vec1,vec2)
% calculate the mutual information of two vectors
% By Hanchuan Peng, April/2002
%

[p12, p1, p2] = PN_estpab(vec1,vec2);
h = PN_estmutualinfo(p12,p1,p2);