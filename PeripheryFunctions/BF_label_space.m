function labels = BF_label_space(k,skip)
% just label with normal labels: everything individually
% only labels every 'skip' objects
% Ben Fulcher 2008/09

if nargin < 2
    skip = 1; % skip not specified: don't skip any
end 
l = length(k);

labels = cell(l,1);
for i = 1:l
    if mod(i,skip) == 0
        labels{i} = k{i};
    else
        labels{i} = '';
    end
end

end