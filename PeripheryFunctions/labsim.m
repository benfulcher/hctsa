function labels = labsim(k,m)
% just label with normal labels: everything individually
% only labels every 'm' objects
% Ben Fulcher 2008/09

if nargin<2; m = 1; end; % m not specified: don't skip any
l = length(k);

labels = cell(l,1);
for i=1:l
    if mod(i,m)==0, labels{i}=k{i};
    else labels{i}='';
    end
end

end