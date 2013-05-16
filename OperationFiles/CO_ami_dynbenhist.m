function out = CO_ami_dynbenhist(y,meth,nbins)
% Look at automutual information using histograms -- vary over number of
% bins to really get a good estimate; return max AMI over bin number range,
% or minimum or median, or....
% Uses hist2 by Nedialko Krouchev (from MATLAB Central)
% [[hist2 for the people by Nedialko Krouchev
%   20 Sep 2006 (Updated 21 Sep 2006)]]
% y should be zscored
% some methods require extra parameter nbins
% Ben Fulcher September 2009

taur = 0:50;
nr = length(taur);
amis = zeros(nr,1);

for i=1:nr;
    amis(i) = CO_ami_benhist(y,taur(i),meth,nbins);
end

% plot(taur,amis,'.-k');


end