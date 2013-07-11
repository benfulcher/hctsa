function out = MS_fnn(y,de,tau,th,kth,justbest,bestp)
% Wrapper for Michael Small's False Nearest Neighbour Code:
% http://small.eie.polyu.edu.hk/matlab/
% Ben Fulcher 19/2/2010

%% INPUTS

% embedding dimension(s), de
if nargin < 2 || isempty(de)
  de = 1:10;
%   disp(['de = ',int2str(de(1)),':',int2str(de(end))]);
end

% Time delay, tau
if nargin < 3 || isempty(tau)
    tau = 1;
%   tau = firstzero(y);
%   disp(['tau = ',int2str(tau)]);
end
if strcmp(tau,'ac')
    tau = CO_fzcac(y); % first zero-crossing of autocorrelation function
elseif strcmp(tau,'mi')
    tau = CO_fmmi(y); % first minimum of automutual information function
end

% A distance threshold for neighbours
if nargin < 4
  th = 5;
%   disp(['th = ',int2str(th)]);
end

% Distance to next points
if nargin < 5
  kth = 1;
%   disp(['kth = ',int2str(kth)]);
end

% (Actually better to use MS_unfolding now -- does a near-identical thing
% to this...)
if nargin < 6 || isempty(justbest)
    justbest = 0;
end
if nargin < 7 || isempty(bestp)
    bestp = 0.05; % first time under 5% of neighest neighbours
end




%% ______________MICHAEL_SMALL'S____CODE__________________
%function [nfnn]=fnn(y,de,tau,th,kth)
%
%determine the number of false nearest neighbours for the time
%series y embedded in dimension de with lag tau. 
%
%for each pair of values (de,tau) the data y is embeded and the
%nearest neighbour to each point (excluding the immediate
%neighbourhood of n points) is determined. If the ratio of the
%distance of the next (kth) points and these points is greater than
%th then they are counted as false nearest neighbours.
%
% default:
% th=5
% kth=1
%
% p(i,j) is the proportion of false nearest neighbours for de(i)
% and tau(j).
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk



p=[];
for t=tau
  px=[];
  for d=de,
    %embed the data
    X = MS_embed(y,d,t); % changed to MS_embed ++BF
    if isnan(X)
        out = NaN;
    end
    [dx,nx]=size(X);

    %find the nearest neighbours of each point
    ind=nearest(X(:,1:(nx-kth)),tau); %whooh hooo!
    

    %distance between each point and its nearest neighbour
    d0=rms(X(:,(1:(nx-kth)))'-X(:,ind)');
    %... and after one time step
    d1=rms(X(:,(kth+1):nx)'-X(:,ind+1)');

    %exclude any coincident points
    d1(d0 == 0)=[];
    d0(d0 == 0)=[];
    
    %calculate the proportion fnn
    ifnn=sum((d1./d0)>th)/length(d0);
    
    %disp -- commented out ++BF
%     disp(['tau = ', int2str(t),', de = ',int2str(d),', nfnn = ',num2str(ifnn*100),'%']);
    

    px=[px ifnn];
  end;
  
  p=[p;px];
    
end;

p=p';

%_____________________________________________________________________

%% Now make output
% Assuming we've set tau, and m is a vector, we should have p (the
% proportion of false neighbours) and de (the corresonding embedding
% dimensions) as vectors


if justbest
    % We just want a scalar to choose the embedding with
    out = de(find(p < bestp,1,'first'));
    if isempty(out), out = de(end)+1; end
    return
else
    % Output all of them
    for i = 1:length(de)
        eval(['out.pfnn_' num2str(de(i)) ' = ' num2str(p(i)) ';']);
    end

    % Output mean
    out.meanpfnn = mean(p);
    % Standard deviation
    out.stdpfnn = std(p);

    % Find when reaches under 20%
    out.firstunder02 = de(find(p<0.2,1,'first'));
    if isempty(out.firstunder02), out.firstunder02 = de(end)+1; end
    % Find when reaches under 10%
    out.firstunder01 = de(find(p<0.1,1,'first'));
    if isempty(out.firstunder01), out.firstunder01 = de(end)+1; end
    % Find first under 5%
    out.firstunder005 = de(find(p<0.05,1,'first'));
    if isempty(out.firstunder005), out.firstunder005 = de(end)+1; end
    % Find first under 2%
    out.firstunder002 = de(find(p<0.02,1,'first'));
    if isempty(out.firstunder002), out.firstunder002 = de(end)+1; end
    % Find first under 1%
    out.firstunder001 = de(find(p<0.01,1,'first'));
    if isempty(out.firstunder001), out.firstunder001 = de(end)+1; end


    % Maximum step-wise change
    out.max1stepchange = max(abs(diff(p)));

end
