function rs = filterbank(s, depth, filterlen)

%tstoolbox/@signal/filterbank
%   Syntax:
%     * filterbank(s, depth, filterlen)
%
%   Filter scalar signal s into 2^textdepth bands of equal bandwith, using
%   maximally flat filters.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,3);

if nargin < 3
	filterlen = 8;
end

nbands = 2^depth;

[h,g]=maxflat(filterlen,filterlen);

c = filterbank(s.core, h, g, 0, depth * ones(1, nbands));

c = core(reshape(data(c), dlens(c)/nbands, nbands));

rs = signal(c, s);	% special constructor calling syntax for working routines
a = getaxis(s, 1);
rs = setaxis(rs, 2, achse(1/unit(a), samplerate(a)/4/nbands, samplerate(a)/2/nbands));
a = setdelta(a, nbands*delta(a));
rs = setaxis(rs, 1, a);
rs = setplothint(rs, 'spectrogram');
rs = addhistory(rs, 'Filterbank');
rs = addcommandlines(rs, 's = filterbank(s', depth, filterlen);




function [h,g,rh,rg]=maxflat(N0,Npi)

%MAXFLAT    Generates maximally flat FIR filters.
%
%           [H,G,RH,RG] = MAXFLAT(N0,NPI) where (N0-1) is the degree
%           of flatness at w=0 and (NPI-1) at w=pi radians. 
%           This function returns half band filters only if N0=NPI.
%           Otherwise, they cannot be passed to WT, WPK.
%
%           H is the analysis lowpass filter, RH the synthesis 
%           lowpass filter, G the analysis highpass filter and
%           RG the synthesis highpass filter.
%
%           For a given order, an increase in NPI results in a wider
%           stopband. 
%
%           References: P.P. Vaidyanathan, "Multirate Systems and
%                       Filter Banks", Prentice-hall, pp. 532-535.

poly=[];                % Calculate trigonometric polynomial
for i=1:N0
        poly=[poly , 2*numcomb(Npi+i-2,i-1)];
end
poly=poly(length(poly):-1:1);
zerospoly=roots(poly);  % Calculate roots

% Transform roots

rootsz=[];

for i=1:length(zerospoly)
    rooty=zerospoly(i);
    rootz1=(1-2*rooty)-2*sqrt(rooty*(rooty-1));
    rootz2=(1-2*rooty)+2*sqrt(rooty*(rooty-1));
    rootsz=[rootsz,rootz1,rootz2];
end     

zeros=rootsz;
N=length(zeros);

% To construct rh for the minimum phase choice, we choose all the zeros 
% inside the unit circle. 

modulus=abs(zeros);

j=1;
for i=1:N
        if modulus(i)<1
                zerosinside(j)=zeros(i);
                j=j+1;
        end
end

An=poly(1);

realzeros=[];
imagzeros=[];
numrealzeros=0;
numimagzeros=0;


Ni=length(zerosinside);

for i=1:(Ni)
        if (imag(zerosinside(i))==0)
                numrealzeros=numrealzeros+1;
                realzeros(numrealzeros)=zerosinside(i);
        else
                numimagzeros=numimagzeros+1;
                imagzeros(numimagzeros)=zerosinside(i); 
                
        end
end

% Construction of rh from its zeros

rh=[1 1];

for i=2:N0
        rh=conv(rh,[1 1]);
end

for i=1:numrealzeros
        rh=conv(rh,[1 -realzeros(i)]);
end

for i=1:2:numimagzeros
        rh=conv(rh,[1 -2*real(imagzeros(i)) abs(imagzeros(i))^2]);
end

% Normalization

rh=sqrt(2)/sum(rh)*rh;

% Calculate h,g,rg

[rh,rg,h,g]=rh2rg(rh);



function [rh,rg,h,g]=rh2rg(rh)

% Calculate rg from rh.

for i=1:length(rh)        
        rg(i) = -(-1)^i*rh(length(rh)-i+1);
end  

% Calculate h and g

h=rh(length(rh):-1:1);
g=rg(length(rg):-1:1);




function y=numcomb(n,k)

% NUMCOMB   combinatorial number.
%
%           NUMCOMB(n,k) calculates the combinatorial number defined
%           as n!/(k!.(n-k)!). 

if n==k,
   y=1;
elseif k==0,
   y=1;
elseif k==1,
   y=n;
else 
   y=fact(n)/(fact(k)*fact(n-k));
end



function y=fact(x)

% FACT   Factorial.
%        FACT(X) is the factorial of the elements in X vector.

for j=1:length(x)
    if x(j)==0,
       y(j)=1;
    else
       y(j)=x(j)*fact(x(j)-1);
    end
end

