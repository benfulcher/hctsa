function out = SY_StatAv(y,n,jo)
% Ben Fulcher 2009
% Ben Fulcher:: minor edits 22/9/2010
% Better to use the 'buffer' function

% jo=='seg' => do it by number of subsegments, n = number of subsegments
% jo=='len' => do it by segment length, n = segment length


if (strcmp(jo,'seg')==1)
    %divide time series into n segments
    M = zeros(n,1);
    p = floor(length(y)/n);% lose the last N mod n data points
    
    for j=1:n
        M(j) = mean(y(p*(j-1)+1:p*j));
    end

elseif (strcmp(jo,'len')==1)
    if length(y) > 2*n
        pn = floor(length(y)/n);
        M = zeros(pn,1);
        for j = 1:pn
            M(j) = mean(y((j-1)*n+1:j*n));
        end
    else
        disp('TS too short for this StatAv')
        out = NaN;
        return
    end
else
    disp('error with StatAv call: either ''seg'' or ''len'' please')
    return
end

s = std(y); % should be 1
sdav = std(M);
out = sdav/s;

end