function out = SY_stdnthder(y,n)
% Idea came from the following source:
% "You can measure the standard deviation of the n-th derivative, if you
% like." -- Vladimir Vassilevsky, DSP and Mixed Signal Design Consultant
% from http://www.mathworks.de/matlabcentral/newsreader/view_thread/136539
% Ben Fulcher, 2010

yd = diff(y,n); % the crudest way of taking a derivative, sure... Could be improved
                % in future
out = std(yd);

end