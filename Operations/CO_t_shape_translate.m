function out = CO_t_shape_translate(y,shape,d,howtomove)
% given a time series, a shape and '~size' d, and a method of moving the shape
% around, returns statistics on how many data points reside inside this
% shape as it is moved around the time series.
% This may be more powerful in an embedding space, but here we do it just
% in the temporal domain (_t_).
% Could also do with a soft boundary, some decaying force function V(r),
% perhaps truncated...?
% Could also do with other shapes -- thick parabolas, squares, ellipses, etc.
% The input time series, y, must be a column vector
% Ben Fulcher, September 2009


%% Check inputs:
if nargin < 2 || isempty(shape)
    shape = 'circle';
end
if nargin < 3 || isempty(d)
    d = 2; % a default distance d = 2
end
if nargin < 4 || isempty(howtomove)
    howtomove = 'pts'; % by default, places shapes on each timepoint
end

N = length(y); % the length of the time series

% y must be a column vector, transpose it if it's a row vector
if size(y,2) > size(y,1)
    y = y';
end

ty = [(1:N)', y]; % has increasing integers as time in the first column

switch howtomove
    case 'pts' % Place shapes on each timepoint (excluding a range at start and end)
        switch shape
            case 'circle' % uses a circle of radius 'd'
                r = d;
                w = floor(r); % only consider a window radius w (these are the
                            %    only points that could possibly be inside)
                rnge = 1+w:N-w;
                NN = length(rnge); % number of admissible points
                np = zeros(NN,1); % number of points
                for i = 1:NN
                    win = ty(rnge(i)-w:rnge(i)+w,:); % create window
                    difwin = win - ones(2*w+1,1)*ty(rnge(i),:);
                    np(i) = sum(sum(difwin.^2,2) <= r^2); % number of points enclosed in shape
                end
                out.max = max(np);
                out.std = std(np);
                
                histnp = zeros(2*w+1,1); % maximum possible hits in circle
                for i = 1:2*w+1
                    histnp(i) = sum(np == i);
                end
                
                [out.npatmode, out.mode] = max(histnp);
                out.npatmode = out.npatmode/NN;
                
                if 2*w + 1 >= 1; out.ones = histnp(1)/NN; end
                if 2*w + 1 >= 2; out.twos = histnp(2)/NN; end
                if 2*w + 1 >= 3; out.threes = histnp(3)/NN; end
                if 2*w + 1 >= 4; out.fours = histnp(4)/NN; end
                if 2*w + 1 >= 5; out.fives = histnp(5)/NN; end
                if 2*w + 1 >= 6; out.sixes = histnp(6)/NN; end
                if 2*w + 1 >= 7; out.sevens = histnp(7)/NN; end
                if 2*w + 1 >= 8; out.eights = histnp(8)/NN; end
                if 2*w + 1 >= 9; out.nines = histnp(9)/NN; end
                if 2*w + 1 >= 10; out.tens = histnp(10)/NN; end
                if 2*w + 1 >= 11; out.elevens = histnp(11)/NN; end
                
                % stationarity in 2,3,4 segments
                % This would be much nicer if I'd have remembered to use
                % the buffer command... :(! -- plus it would have been better (equal segment sizes...)
                div2 = round(linspace(1,NN,3));
                div3 = round(linspace(1,NN,4));
                div4 = round(linspace(1,NN,5));
                
                out.statav2_m = std([mean(np(div2(1):div2(2))), mean(np(div2(2)+1:div2(3)))])/std(np);
                out.statav2_s = std([std(np(div2(1):div2(2))), std(np(div2(2)+1:div2(3)))])/std(np);
                
                out.statav3_m = std([mean(np(div3(1):div3(2))), mean(np(div3(2)+1:div3(3))), ...
                                    mean(np(div3(3)+1:div3(4)))])/std(np);
                out.statav3_s = std([std(np(div3(1):div3(2))), std(np(div3(2)+1:div3(3))), ...
                                    std(np(div3(3)+1:div3(4)))])/std(np);
                
                out.statav4_m = std([mean(np(div4(1):div4(2))), mean(np(div4(2)+1:div4(3))), ...
                                  mean(np(div4(3)+1:div4(4))), mean(np(div4(4)+1:div4(5)))])/std(np);
                out.statav4_s = std([std(np(div4(1):div4(2))), std(np(div4(2)+1:div4(3))), ...
                                  std(np(div4(3)+1:div4(4))), std(np(div4(4)+1:div4(5)))])/std(np);
        otherwise
            error('Unknwon shape ''%s''',shape)
        end
    otherwise
        error('Unknwon setting for ''howtomove'' input: ''%s''',howtomove)
end

% plot(np)



end