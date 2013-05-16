function out = IN_cgent(y,ng,howtocg)
% Determines the entropy given a course-graining
% Ben Fulcher, September 2009

% Doesn't really make sense since quantiles gives each letter in the
% alphabet an equal probability... Entropy is thus signal independent.

% plot(y)
% % Course-grain the time series
% yth = SUB_coursegrain(y,ng,howtocg);
% hold on; plot(zscore(yth),':r');
% N=length(yth);
% 
% % determine probability of each letter in alphabet
% ps=zeros(ng,1);
% for i=1:ng
%     ps(i)=length(find(yth==i));
% end
% ps
% ps=ps/N;
% 
% out=-sum(ps.*log(ps));

end