function entropy = EN_pdfent(y)

tmp = sort(y);
tmp = diff(tmp);
if size(tmp,1)<size(tmp,2)
    tmp=tmp';
end
pdf = diff([0; find(tmp>0); length(y)])/length(y); % empirical probability distribution
entropy = -sum(pdf.*log2(pdf));

end