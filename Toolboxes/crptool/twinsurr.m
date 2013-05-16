function [y nTw] = twinsurr(varargin)
% TWINSURR   Creates twin surrogates for statistical tests.
%    Y=TWINSURR(X) creates twin surrogates Y based on the vector X
%    using recurrences. The matrix Y contains 100 columns of
%    100 twin surrogates. If X is a PxQ matrix, the resulting
%    surrogate matrix is Px100xQ.
%
%    Y=TWINSURR(X,M,T,E,...) creates twin surrogates using
%    embedding dimension M, delay T, recurrence threshold E. The
%    input arguments are similar to those of the command CRP.
%
%    Y=TWINSURR(X,M,T,E,...,N) creates N surrogates (default is 100).
%
%    [Y,NTWINS]=TWINSURR(...) where NTWINS is the total number
%    of twins in the RP.
%
%    Example: x = rand(3,1);
%             a = [.8 .3 -.25 .9]';
%             for i = 4:1000,
%                x(i) = sum(a(1:3) .* x(i-1:-1:i-3)) + a(end) * randn;
%             end
%             xs = twinsurr(x,1,1,.1,'euc',10);
%
%    See also CRP, RECONS.
%
%    References: 
%    Thiel, M., Romano, M. C., Kurths, J., Rolfs, M., Kliegl, R.: 
%    Twin Surrogates to Test for Complex Synchronisation, 
%    Europhys. Lett., 75, 2006.

% Copyright (c) 2008 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2008/07/02 12:00:48 $
% $Revision: 5.2 $
%
% $Log: twinsurr.m,v $
% Revision 5.3  2008/09/24 12:00:48  marwan
% serious bug fix, zero columns were also considered as twins!
%
% Revision 5.2  2008/07/02 12:00:48  marwan
% silent ability added, minor bug fixes
%
% Revision 5.1  2008/07/01 13:09:27  marwan
% initial import
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% programme properties

global errcode props

init_properties
nsur_init = 100;
sil = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check the input

error(nargchk(1,7,nargin));
if nargout>2, error('Too many output arguments'), end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read the input

% transform any int to double
intclasses = {'uint8';'uint16';'uint32';'uint64';'int8';'int16';'int32';'int64'};
flagClass = [];
for i = 1:length(intclasses)
    i_int=find(cellfun('isclass',varargin,intclasses{i}));
    if ~isempty(i_int)
        for j = 1:length(i_int)
            varargin{i_int(j)} = double(varargin{i_int(j)});
        end
        flagClass = [flagClass; i_int(:)];
    end
end
if ~isempty(flagClass)
    disp(['Warning: Input arguments at position [',num2str(flagClass'),'] contain integer values']);
    disp(['(now converted to double).'])
end


varargin{8}=[];
i_double=find(cellfun('isclass',varargin,'double'));
i_char=find(cellfun('isclass',varargin,'char'));
check_meth={'ma','eu','mi','nr','rr','fa','in','om','op','di'}; 	% maxnorm, euclidean, nrmnorm,  fan, distance
check_norm={'non','nor'};                        % nonormalize, normalize
check_sil={'ve','si'};                         % verbose, silent

if isnumeric(varargin{1}) 		% read commandline input
   % check the text input parameters for method, gui 
    temp_meth=0;
    temp_norm=0;
    temp_sil=0;
    if ~isempty(i_char)
         for i=1:length(i_char), 
            varargin{i_char(i)}(4)='0';
            temp_norm=temp_norm+strcmpi(varargin{i_char(i)}(1:3),check_norm'); 
            temp_meth=temp_meth+strcmpi(varargin{i_char(i)}(1:2),check_meth'); 
            temp_sil=temp_sil+strcmpi(varargin{i_char(i)}(1:2),check_sil'); 
         end
         method_n=min(find(temp_meth));
         nonorm=min(find(temp_norm))-1;
         sil=min(find(temp_sil))-1;
         for i=1:length(i_char); temp2(i,:)=varargin{i_char(i)}(1:3); end

         if isempty(sil), sil=0; end
         if isempty(nonorm), nonorm=1; end
         if nonorm>1, nonorm=1; end
         if isempty(method_n), method_n=1; end
         if method_n>length(check_meth), method0=length(check_meth); end
         method=check_meth{method_n};
         norm_str = check_norm{nonorm+1};
    else
         method = 'max'; norm_str = 'non';
    end

    % get the parameters for creating RP
    if max(size(varargin{1}))<=3
        disp('Error using ==> twinsurr')
        disp('To less values in data X.')
        return
    end
    x=double(varargin{1});

    if ~isempty(varargin{i_double(2)}), m=varargin{i_double(2)}(1); else m=1; end
    if ~isempty(varargin{i_double(3)}), t=varargin{i_double(3)}(1); else t=1; end
    if ~isempty(varargin{i_double(4)}), e=varargin{i_double(4)}(1); else e=.1; end
    if ~isempty(varargin{i_double(5)}), nsur=varargin{i_double(5)}(1); else nsur=nsur_init; end
else
    disp('Error using ==> twinsurr')
    disp('No valid arguments.')
    return
end


if size(x,1) < size(x,2), x = x'; end
m0 = size(x,2);
N = length(x); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% make surrogates

X = crp(x,m,t,e,method,norm_str,'sil');
X = double(X) .* (1-eye(size(X)));
NX = length(X);
%% find twins
if ~sil, h = waitbar(0,'Searching Twins'); end
for i = 1:NX, if ~sil, waitbar(i/NX); end
    A = repmat(X(:,i),1,NX);
%    S{i} = find(all(X == A & X == 1));
    S{i} = [i find(all(X == A) & any(X(:,i)))];
    nTwins(i) = numel(S{i});
end
if ~sil, delete(h); end

nTwins = sum(nTwins) - NX;

if nTwins < 1
    warning(['No twins found!',10,'The derived surrogates are identical with the original data.']')
end


for k = 1:nsur
    %% chose randomly first point
    i = ceil(NX * rand);
    xs = x(i,:);

    %% jump randomly between twins until length of surrogate is reached
    while length(xs) < N
        j = S{i};
        j_i = ceil(length(j) * rand);
        i = j(j_i);
        i = i + 1;
        if i > NX
            i = ceil(NX * rand);
            continue 
        end
        xs = [xs; x(i,:)];
    end
    y(:,k,1:m0) = reshape(xs(:),N,1,m0);
end

if nargout == 2;
    nTw = nTwins;
end
