function cout=filterbank(cin,h,g,order,basis)

%tstoolbox/@core/filterbank
%   Syntax:
%     * filterbank(cin,H,G,ORDER,BASIS)
%
%   Input Arguments:
%     * H - lowpass filter
%     * G - highpass filter
%     * ORDER - indicates the type of tree:
%          + 0 - band sorting according to the filter bank
%          + 1 - band sorting according to the frequency decomposition
%     * BASIS - desired subband decomposition
%
%   calculates the Wavelet Packet Transform of cin. It can be obtained
%   using a selection algorithm function. It may be switched from one
%   format to another using CHFORMAT. The different bands are sorted
%   according to ORDER and BASIS. If BASIS is omitted, the output is a
%   matrix with the coefficients obtained from all the wavelet packet
%   basis in the library. Each column in the matrix represents the outputs
%   for a level in the tree. The first column is the original signal. If
%   the length of X is not a power of 2, the columns are zero padded to
%   fit the different lengths. Run the script 'BASIS' for help on the
%   basis format.
%   See also: IWPK, CHFORMAT, PRUNEADD, PRUNENON, GROWADD, GROWNON.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

%       Auxiliary input argument for recursion:
%       flag: if set to 1, indicates that no basis is specified
%          and the output will be a matrix
%       flag2: if set to 1 the function was called with the detail signal,
%          then the coefficients must be rearranged

cout = core(wpk(data(cin), h,g,order,basis));


function w = wpk(x, h,g,order,basis,flag,flag2)

if ((nargin==4)|(nargin==5))
        if nargin==4
                L=floor(log2(length(x)));       % initialization when
                basis=ones(2^L,1)*(L-1);        % no basis was specified
                flag=1;
                flag2=0;
        else                                    % basis was specified 
                if (~basis)
                        w=x;
                        return
                end
                basis=basis-1;
                flag=0;
                flag2=0;
        end
        if size(x,2)~=1                 % shapes the signal as column
                x=x(:);                 % vector if necessary
                fil=1;
        else
                fil=0;
        end
end

wx=wt(x,h,g,1);                         % perform one analysis level
                                        % into the analysis tree

N=length(wx);                           % separate approximation and detail
a=wx(1:N/2);                            % at the analysis output
d=wx(N/2+1:N);

if (flag2*order)                        % if the recursion is going to be
        s1=d;s2=a;                      % performed on the detail signal,
else                                    % and it is frequency sorted, the
        s1=a;s2=d;                      % outputs must be swapped
end

if all(basis==0)                           % two ending nodes achieved
        w=[s1;s2];                         
        if (nargin==5)&(fil), w=w'; end
        return
end;

tope=1;                                 % finds the point where the
suma=2^(-basis(1));                     % basis vector must be divided
i=1;                                    % in order to call the recursion
while (suma<tope)
        i=i+1;
        suma=suma+2^(-basis(i));
end

if all(basis(1:i)~=0)                                   % a level with an ending node
        wa=wpk(s1,h,g,order,(basis(1:i)-1),flag,0);     % but the other node continues
else wa=s1;
end

if all(basis(i+1:length(basis))~=0)
        wd=wpk(s2,h,g,order,(basis(i+1:length(basis))-1),flag,1);   % complementary case
else wd=s2;
end

if flag                                 % output is a matrix
        lwd=length(wd);
        lwx=length(wx);
        nz=2*lwd-lwx;                   % find the number of zeros
                                        % to be inserted for padding

        w=[ [s1;zeros(nz/2,1);s2;zeros(nz/2,1)] , [wa;wd] ];
else
        w=[wa;wd];
end

if nargin==4
        N=size(w,1);
        w=[ [x;zeros(N-length(x),1)] w];
elseif (nargin==5)&(fil), w=w'; end



function y=wt(x,h,g,k,del1,del2)

% WT   Discrete Wavelet Transform.
% 
%      WT(X,H,G,K) calculates the wavelet transform of vector X. 
%      If X is a matrix (2D), WT will calculate the one dimensional 
%      wavelet transform of each row vector. The second argument H 
%      is the lowpass filter and the third argument G the highpass
%      filter.
%
%      The output vector contains the coefficients of the DWT ordered 
%      from the low pass residue at scale K to the coefficients
%      at the lowest scale, as the following example ilustrates:
%
%      Output vector (k=3):
%
%      [------|------|------------|------------------------]
%         |       |        |                |
%         |       |        |                `-> 1st scale coefficients 
%         |       |        `-----------> 2nd scale coefficients
%         |       `--------------------> 3rd scale coefficients
%         `----------------> Low pass residue  at 3rd scale 
%
%       
%      If X is a matrix, the result will be another matrix with 
%      the same number of rows, holding each one its respective 
%      transformation.
%
%      WT (X,H,G,K,DEL1,DEL2) calculates the wavelet transform of 
%      vector X, but also allows the users to change the alignment
%      of the outputs with respect to the input signal. This effect
%      is achieved by setting to DEL1 and DEL2 the delays of H and
%      G respectively. The default values of DEL1 and DEL2 are 
%      calculated using the function WTCENTER. 

% -----------------------------------
%    CHECK PARAMETERS AND OPTIONS
% -----------------------------------

h=h(:)';        % Arrange the filters so that they are row vectors.
g=g(:)';

if length(x)<2^k 
        disp('The scale is too high. The maximum for the signal is:')
        floor(log2(length(x)))
        return
end

[liy,lix]=size(x);

if lix==1               % And arrange the input vector to a row if 
        x=x';           % it's not a matrix. 
        trasp=1;        % (and take note of it)
        [liy,lix]=size(x);
else
        trasp=0;
end


%--------------------------
%    DELAY CALCULATION 
%--------------------------

% Calculate delays as the C.O.E. of the filters
dlp=wtcenter(h);
dhp=wtcenter(g);

if rem(dhp-dlp,2)~=0            % difference between them.
        dhp=dhp+1;              % must be even
end;

if nargin==6,                   % Other experimental filter delays
        dlp=del1;               % can be forced from the arguments
        dhp=del2;
end;

%------------------------------
%    WRAPPAROUND CALCULATION 
%------------------------------
llp=length(h);                  % Length of the lowpass filter
lhp=length(g);                  % Length of the highpass filter.

L=max([lhp,llp,dlp,dhp]);       % The number of samples for the
                                % wrapparound. Thus, we should need to 
                                % move along any L samples to get the
                                % output wavelet vector phase equal to
                                % original input phase.


%------------------------------
%     START THE ALGORITHM 
%------------------------------

for it=1:liy,           % For every row of the input matrix...
                        % (this makes one wavelet transform
                        % for each of the rows of the input matrix)
        tm=[];
        t=x(it,:);                      % Copy the vector to transform.

        for i=1:k                       % For every scale (iteration)...
                lx=length(t);
                if rem(lx,2)~=0         % Check that the number of samples
                        t=[t,0];        % will be even (because of decimation).
                        lx=lx+1;
                end
                tp=t;                   % Build wrapparound. The input signal
                pl=length(tp);          % can be smaller than L, so it can
                while L>pl              % be necessary to repeat it several
                        tp=[tp,t];      % times
                        pl=length(tp);
                end

                t=[tp(pl-L+1:pl),t,tp(1:L)];    % Add the wrapparound.

                yl=conv(t,h);           % Then do lowpass filtering ...
                yh=conv(t,g);           % ... and highpass filtering.

                yl=yl((dlp+1+L):2:(dlp+L+lx));    % Decimate the outputs
                yh=yh((dhp+1+L):2:(dhp+L+lx));    % and leave out wrapparound

                tm=[yh,tm];             % Put the resulting wavelet step
                                        % on its place into the wavelet 
                                        % vector...
                t=yl;                   % ... and set the next iteration.
        end

        y(it,:)=[t,tm];                 % Wavelet vector (1 row vector)


end                             % End of the "rows" loop.

%------------------------------
%    END OF THE ALGORITHM 
%------------------------------

if trasp==1                     % If the input data was a column vector
        y=y';                   % then transpose it.
end



function d=wtcenter(x,op);

%  WTCENTER Calculates the delay of filters for alignment.
%
%           WTCENTER (X) calculates the integer aproximation
%           of delay for filter X using the method set with
%           the WTMETHOD function, for alignment operations
%           in Wavelet transforms.
%
%           For a non integer value, use the CENTER function.
%
%           See also: WTMETHOD, CENTER, WT
%

global WTCENTERMETHOD

if size(WTCENTERMETHOD)==[0,0]
        WTCENTERMETHOD=0;
end

if WTCENTERMETHOD>3 | WTCENTERMETHOD<0
        WTCENTERMETHOD=0
end
        
d=floor(center(x,WTCENTERMETHOD));

% (Another long function !!!)


function d=center(x,op);

%  CENTER  Delay calculation for Wavelet transform alignment.
%
%          CENTER (X, OP) calculates the delay for filter in X 
%          according to the alignment method indicated in OP. 
%          This delay is used by Wavelet transform functions.
%          The value of OP can be:
%              0 : First Absolute Maxima Location
%              1 : Zero delay in analysis (Full for synthesis).
%              2 : Mass center (sum(m*d)/sum(m))
%              3 : Energy center (sum(m^2 *d)/sum(m^2))
%
%          If no output argument is given, then the vector X will
%          be plotted in the current figure, and a color line will be 
%          marking the result.(red: OP=0; green: OP=1; cyan: OP=2; 
%          blue: OP=4)

lx=length(x);
l=1:lx;

if op==1
        d=0;
else
        if op==2
                xx=abs(x(:)');
                L=l;
        end
        if op==3
                xx=x(:)'.^2;
                L=l;
        end
        if op==0
                [mx,d]=max(abs(x));
        else 
                
                d=sum(xx.*L)/sum(xx);
        end
end


