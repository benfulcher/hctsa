function ref = randref(low, high, N)

% Create roughly N random integers between between
% low and high. Each value is unique, the values
% are monotonically increasing.
% These integer values can be used as reference indices
% for statistical purposes
% -1 for N means : return low:high
%
% length(randref(low, high, N)) will not be exactly equal to N !!!

narginchk(3,3);

if (high<=low)
	error('Cannot create random indices when upper limit is smaller than lower limit');
end

if N < 1 | N > high-low
    ref = low:high;
    return
end

m_target = (high-low) / N;

if m_target < 10
    m = (0.99+log10(m_target)) *  m_target;
    N  = ceil(1.5 * N);     % calculate more indices than necesary
else
    m = 1.99 * m_target;
    N = ceil(1.2*N);        % calculate more indices than necesary
end

ref = cumsum([low ; ceil(rand(N,1)*m)]);    % compute indices
ref = ref(find(ref <= high));               % remove out of range indices 
