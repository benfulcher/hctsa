function out = SY_stdnthder(y,n)
% Idea came from the following source:
% "You can measure the standard deviation of the n-th derivative, if you
% like." -- Vladimir Vassilevsky, DSP and Mixed Signal Design Consultant
% from http://www.mathworks.de/matlabcentral/newsreader/view_thread/136539
% Ben Fulcher, 2010

if nargin < 2 || isempty(n)
    n = 2;
end

yd = diff(y,n); % crude method of taking a derivative that could be improved
                % upon in future
out = std(yd);

end