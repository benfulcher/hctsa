function out = SY_StatAv(y,n,whattype)
% Ben Fulcher 2009
% Ben Fulcher:: minor edits 22/9/2010
% Better to use the 'buffer' function

% whattype == 'seg' => do it by number of subsegments, n = number of subsegments
% whattype == 'len' => do it by segment length, n = segment length

if nargin < 2
    n = 5; % use 5 segments
end
if nargin < 3
    whattype = 'seg'; % divide into n segments by default
end

N = length(y); % time-series length

switch whattype
case 'seg'
    % divide time series into n segments
    M = zeros(n,1);
    p = floor(N/n);% lose the last N mod n data points
    
    for j = 1:n
        M(j) = mean(y(p*(j-1)+1:p*j));
    end

case 'len'
    if N > 2*n
        pn = floor(N/n);
        M = zeros(pn,1);
        for j = 1:pn
            M(j) = mean(y((j-1)*n+1:j*n));
        end
    else
        fprintf(1,'This time series (N = %u) is too short for StatAv(%s,%u)\n',N,whattype,n)
        out = NaN; return
    end
otherwise
    error('Error evaluating StatAv of type ''%s'', please select either ''seg'' or ''len''',whattype)
end

s = std(y); % should be 1 (for a z-scored time-series input)
sdav = std(M);
out = sdav/s;

end