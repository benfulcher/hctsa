% ZG_rsum(X)
% row sum
% 
% Machine Learning Toolbox
% Version 1.0  01-Apr-96
% Copyright (c) by Zoubin Ghahramani
% http://mlg.eng.cam.ac.uk/zoubin/software.html
%
% ------------------------------------------------------------------------------
% The MIT License (MIT)
% 
% Copyright (c) 1996, Zoubin Ghahramani
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
% ------------------------------------------------------------------------------

function Z = ZG_rsum(X)

Z = zeros(size(X(:,1)));

for i = 1:length(X(1,:))
    Z = Z + X(:,i);
end

end