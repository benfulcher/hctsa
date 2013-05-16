function out = SY_stdnthder(y,n)
% "You can measure the standard deviation of the n-th derivative, if you
% like." -- Vladimir Vassilevsky, DSP and Mixed Signal Design Consultant
% from http://www.mathworks.de/matlabcentral/newsreader/view_thread/136539

yd = diff(y,n); % the coarsest way of taking a derivative, sure... Could be improved
              % in future
out = std(yd);


end