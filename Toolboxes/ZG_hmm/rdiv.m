% function Z=rdiv(X,Y)
%
% row division: Z = X / Y row-wise
% Y must have one column 

function Z=rdiv(X,Y)

if(length(X(:,1)) ~= length(Y(:,1)) | length(Y(1,:)) ~=1)
  disp('Error in RDIV');
  return;
end

Z=zeros(size(X));

for i=1:length(X(1,:))
  Z(:,i)=X(:,i)./Y;
end
