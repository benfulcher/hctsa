function h = PN_entropy(vec1)
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
% h = entropy(vec1)
% calculate the entropy of a variable vec1
%
% demo: 
%  a=[1 2 1 2 1]';b=[2 1 2 1 1]';
%  fprintf('entropy(a) = %d\n',entropy(a));
%
% the same as entropycond(vec1)
%
% By Hanchuan Peng, April/2002
%

%% BEN -- REMOVED THIS ERROR CHECKING, I'D PREFER THE ERROR
% if nargin < 1,
% 
%   disp('Usage: h = entropy(vec1).');
%   h = -1;
% 
% else,

  [p1] = PN_estpa(vec1);
  h = PN_estentropy(p1);

% end;


