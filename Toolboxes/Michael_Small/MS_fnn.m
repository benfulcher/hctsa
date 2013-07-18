function [p]=fnn(y,de,tau,th,kth);
  
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

if nargin<5,
  kth=1;
  disp(['th = ',int2str(kth)]);
end;

if nargin<4,
  th=5;
  disp(['th = ',int2str(th)]);
end;
if nargin<3,
  tau=MS_firstzero(y);
  disp(['tau = ',int2str(tau)]);
end;

if nargin<2,
  de=[1:10];
  disp(['de = ',int2str(de(1)),':',int2str(de(end))]);
end;

p=[];
for t=tau,
  px=[];
  for d=de,
    %embed the data
    X=MS_embed(y,d,t);
    [dx,nx]=size(X);

    %find the nearest neighbours of each point
    ind=MS_nearest(X(:,1:(nx-kth)),tau); %whooh hooo!
    

    %distance between each point and its nearest neighbour
    d0=MS_rms(X(:,(1:(nx-kth)))'-X(:,ind)');
    %... and after one time step
    d1=MS_rms(X(:,(kth+1):nx)'-X(:,ind+1)');

    %exclude any coincident points
    d1(d0==0)=[];
    d0(d0==0)=[];
    
    %calculate the proportion fnn
    ifnn=sum((d1./d0)>th)/length(d0);
    
    %disp
    % disp(['tau = ', int2str(t),', de = ',int2str(d),', nfnn = ',num2str(ifnn*100),'%']);
    
    px=[px ifnn];
  end;
  
  p=[p;px];
    
end;

p=p';
