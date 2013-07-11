function [chsq,prob,nu]=chisq(b1,b2,nu);

% function [chsq,prob,nu]=chisq(b1,b2,nu);
%
% given the two matrices b1 and b2 of binned data b1(:,1) are the bin edges 
% and b1(:,2) are the number of observations in each bin. This thingo calculates
% the chi-square statistic (chsq);
%           sum(((b1-b2)^2)/(b1b2))
% and the probability (prob) of observing the above value of chi^2 given the
% null hypothesis that the distributions ARE the same. These formulas are from
% "Numerical Recipes in C" pgs. 448-490 and 171-177.
% Calculating this probability requires the number of degrees of freedom, df.
% This function takes 
%         df=(# of bins for which either b1(:,2) or b2(:,2) is positive)-1 
% if no alternative is given.
%
% If the binning is not the same in both cases (i.e. b1(:,1)~=b2(:,1) ) then
% this function assumes b2(:,2)=f(b2(:,1)) for some function f and performs a 
% linear interpolation so that the new distribution bb2 for b2 is 
% b2(:,1)=b1(:,1) and b2(:,2)~f(b1(:,1)). Of course some care needs to be taken
% to evaluate both functions only over the interval which is common to both
% binnned distributions.
%
% This algorithm works fine (ie, is intended to work) with the output from the
% binner module in dimfit.
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk

if nargin < 3,
   nu=[];
end;

%best to just check that both b1(:,1) and b2(:,1) are monotonic. If not we
%will be rather callous and just get rid of the offending bits.
Ind=find(diff(b1(:,1))<=0);
if ~isempty(Ind)
	disp('WARNING binning is stuffed - destroying bins... results may be crap.'); 
	b1=b1(find(diff(b1(:,1))>0),:);
end;
Ind=find(diff(b2(:,1))<=0);
if ~isempty(Ind)
	disp('WARNING: binning is stuffed - destroying bins... results may be crap.'); 
	b2=b2(find(diff(b2(:,1))>0),:);
end;

%well, now we have to make them the same so we better find the common 
%interval.
me=max(min(b1(:,1)),min(b2(:,1)));
Me=min(max(b1(:,1)),max(b2(:,1)));
%now check that there is some overlap
if me>=Me 
	chsq=nan; prob=nan; nu=nan;
	disp('ERROR: intervals do not overlap in chisq');
	return;
end;
%and then look for the points closest to (but not outside of) me and Me in b1.
mi=min(find(b1(:,1)>=me)); me=b1(mi,1);
Mi=max(find(b1(:,1)<=Me)); Me=b1(Mi,1);
ei=mi:1:Mi; bb1=b1(ei,:);

%now we gotta make approximations to b2(:,2).
bb2(:,1)=bb1(:,1);
bb2(:,2)=interp1(b2(:,1),b2(:,2),bb1(:,1));

%we should just check that there are no bins with negative occupancies.
if ~isempty(find(bb1(:,2)<0 | bb2(:,2)<0))
	disp('You moron - some of the bins have less than zero elements.');
	disp('Any result you get from this algorithm will be crap.');
end;

%plot(b1(:,1),b1(:,2)); hold on;
%plot(bb1(:,1),bb1(:,2),'y+');
%plot(b2(:,1),b1(:,2),'g');
%plot(bb2(:,1),bb2(:,2),'g+');

%and I gues we are now only really interested in the bin occupancies.
bb1=bb1(:,2); bb2=bb2(:,2);

%
%
%now for the actual calculation
%

Ind=find((bb1+bb2)>0);
if isempty(nu),
  nu=length(Ind)-1;
end;
chsq=sum(((bb1(Ind)-bb2(Ind)).^2)./(bb1(Ind)+bb2(Ind)));
prob=1-gammainc(chsq/2,nu/2);
