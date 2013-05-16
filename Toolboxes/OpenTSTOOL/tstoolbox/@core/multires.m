function cout=multires(cin,h,rh,g,rg,sc);

%tstoolbox/@core/multires
%   Syntax:
%     * multires(cin,h,rh,g,rg,sc)
%
%   Input Arguments:
%     * cin - core object
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


t=data(cin);
h=h(:);
rh=rh(:);
g=g(:);
rg=rg(:);

lx=length(t);
lf=length(h)+length(rh);
d(1)=floor((lf)/2)-1;

if sc>1
	for i=2:sc
		d(i)=d(i-1)+d(1)*2^(i-1);
	end
end

lt=lx;

y=zeros(lx,sc+1);

tt=t;

for i=1:sc,
        tg=conv(g,tt);
        tg=tg(1:2:length(tg));
        tt=conv(h,tt);
        tt=tt(1:2:length(tt));
        for j=i:-1:1
                if j==i
                        tm=[tg' ; zeros(1, length(tg))];
                        tm=tm(:);
                        tm=conv(rg,tm);
                else
                        tm=[tm' ; zeros(1, length(tm))];
                        tm=tm(:);
                        tm=conv(rh,tm);
                end
        end
        y(:,i)=tm(1+d(i):lx+d(i));              
end

for j=1:sc
        tt=[tt' ;zeros(1, length(tt))];
        tt=tt(:);
        tt=conv(rh,tt);
end

y(:,sc+1)=tt(1+d(sc):lx+d(sc));

cout = core(y);

