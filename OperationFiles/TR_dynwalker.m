function out=TR_dynwalker(y,wmeth,prange)
% looks at statistics from TR_walker as a function of some parameter
% prange should be a vector relevant to the method wmeth

fouts=cell(length(prange),1);
for i=1:length(prange)
    fouts{i}=cell2mat(struct2cell(TR_walker(y,wmeth,prange(i))));
end
% fouts is a cell of vectors
fouts=cell2mat(fouts');
% fouts is a matrix containing the variation of tests across the entire set
% in TR_walker.

out=fouts';
% plot(fouts')

end