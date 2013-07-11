function out = SY_compdtn(y,nseg,eachorpar)
% compares the distribution in n consecutive partitions of the signal,
% returning the sum of differences of the kernel smoothed distributions
% eachorbet is either 'each': compares each subdistribution to each other
% subdistribtuion; or 'par': compares each subdistribtuion to the parent
% (full signal) distribution
% Ben Fulcher August 2009

N = length(y);
lseg = floor(N/nseg);
% z=zeros(lseg,1);
dns = zeros(200,nseg);
r = linspace(min(y),max(y),200); % make range of ksdensity uniform across all subsegments
for i = 1:nseg
    dns(:,i) = ksdensity(y((i-1)*lseg+1:i*lseg),r,'function','pdf');
end
% plot(dns)
% disp('done')

switch eachorpar
    case 'par'
        pardn = ksdensity(y,r,'function','pdf');
        divs = zeros(nseg,1);
        for i = 1:nseg;
            divs(i) = sum(abs(dns(:,i)-pardn')); % each is just divergence to parent
        end
%         hold on;plot(pardn,'r','LineWidth',2);hold off
        % return same statistics as for the 'each' case
        out.meandiv = mean(divs);
        out.mediandiv = median(divs);
        out.mindiv = min(divs);
        out.maxdiv = max(divs);
        out.stddiv = std(divs);
        
    case 'each'
        if nseg == 2 % output is just an integer: only two distributions to compare
            out = sum(abs(dns(:,1)-dns(:,2)));
        else % need to compare a number of different distributions against each other
            diffmat = zeros(nseg); % store pairwise differences
            for i = 1:nseg
                for j = 1:nseg
                    if j>i
                        diffmat(i,j) = sum(abs(dns(:,i)-dns(:,j))); % store sum of absolute differences
                    end
                end
            end
            
            divs = diffmat(diffmat>0); % a vector of all divergences
            if isempty(divs);
                disp('What!? No divergences?! Returning NaNs...');
                out = NaN; return
            end
            out.meandiv = mean(divs);
            out.mediandiv = median(divs);
            out.mindiv = min(divs);
            out.maxdiv = max(divs);
            out.stddiv = std(divs);
        end
        
end

end