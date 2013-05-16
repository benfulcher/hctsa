function varargout = rrspec(varargin)
% RRSPEC   Tau-recurrence rate spectrum.
%    RRSPEC(X,M,T,E,W,FS,...) calculates the tau-recurrence rate
%    spectrum based on a recurrence plot using embedding dimension
%    M, embedding delay T, recurrence threshold E, maximal lag
%    for tau-recurrence W, and sampling frequency FS. The
%    input arguments are similar to those of the command TAUCRP.
%
%    P = RRSPEC(...) returns the tau-recurrence rate spectrum
%    in vector P.
%
%    [P F] = RTSPEC(...) returns the tau-recurrence rate spectrum
%    in vector P and the vector of corresponding frequencies F.
%
%    Example: fs = 22;
%             x = sin(2*pi * [0:1/fs:44]);
%             rrspec(x,2,1,.1,[],fs)
%
%    See also TAUCRP, RTSPEC.
%    

% Copyright (c) 2008 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2008/07/02 12:01:13 $
% $Revision: 5.1 $
%
% $Log: rrspec.m,v $
% Revision 5.1  2008/07/02 12:01:13  marwan
% initial import
%
% Revision 5.1  2008/07/01 13:09:27  marwan
% initial import
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% programme properties

global errcode props

init_properties
fs_init = 100;
w_init = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check the input

error(nargchk(1,8,nargin));
if nargout>2, error('Too many output arguments'), end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read the input

% transform any int to double
intclasses = {'uint8';'uint16';'uint32';'uint64';'int8';'int16';'int32';'int64';'logical'};
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


varargin{9}=[];
i_double=find(cellfun('isclass',varargin,'double'));
i_char=find(cellfun('isclass',varargin,'char'));
check_meth={'ma','eu','mi','nr','rr','fa','in','om','op','di'}; 	% maxnorm, euclidean, nrmnorm,  fan, distance
check_norm={'non','nor'};                        % nonormalize, normalize

if isnumeric(varargin{1}) 		% read commandline input
   % check the text input parameters for method, gui 
    temp_meth=0;
    temp_norm=0;
    if ~isempty(i_char)
         for i=1:length(i_char), 
            varargin{i_char(i)}(4)='0';
            temp_norm=temp_norm+strcmpi(varargin{i_char(i)}(1:3),check_norm'); 
            temp_meth=temp_meth+strcmpi(varargin{i_char(i)}(1:2),check_meth'); 
         end
         method_n=min(find(temp_meth));
         nonorm=min(find(temp_norm))-1;
         for i=1:length(i_char); temp2(i,:)=varargin{i_char(i)}(1:3); end

         if isempty(nonorm), nonorm=1; end
         if nonorm>1, nonorm=1; end
         if isempty(method_n), method_n=1; end
         if method_n>length(check_meth), method0=length(check_meth); end
         method=check_meth{method_n};
         norm_str = check_norm{nonorm+1};
    else
         method = 'max'; norm_str = 'nor';
    end

    % get the parameters for creating RP
    if max(size(varargin{1}))<=3
        disp('Error using ==> rrspec')
        disp('To less values in data X.')
        return
    end
    x=double(varargin{1});

    if ~isempty(varargin{i_double(2)}), m=varargin{i_double(2)}(1); else m=1; end
    if ~isempty(varargin{i_double(3)}), t=varargin{i_double(3)}(1); else t=1; end
    if ~isempty(varargin{i_double(4)}), e=varargin{i_double(4)}(1); else e=.1; end
    if ~isempty(varargin{i_double(5)}), w=varargin{i_double(5)}(1); else w=w_init; end
    if ~isempty(varargin{i_double(6)}), fs=varargin{i_double(6)}(1); else fs=fs_init; end
else
    disp('Error using ==> rrspec')
    disp('No valid arguments.')
    return
end


if size(x,1) < size(x,2), x = x'; end
N = length(x)-(m-1)*t;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% make recurrence spectrum

%% RP
X = taucrp(x,m,t,e,N,method,norm_str,'sil');

%% tau-RR
scaling = [1:N N+1 N:-1:1];
tauRR = sum(X') ./ scaling;

%% Fourier transformation
fft_RR = fft(tauRR-mean(tauRR));
%fft_RR = fhtseq(tauRR-mean(tauRR)); % Walsh spectrum >>> experimental

%% spectrum
P = fft_RR .* conj(fft_RR);
%P = fft_RR.^2;
f = fs*linspace(0,.5,N);
P = P(1:N);

if nargout == 1
    varargout{1} = P(:);
elseif nargout == 2
    varargout{1} = P(:);
    varargout{2} = f(:);
else
    %% Plot the spectrum
    semilogy(f,P(1:N)+1)
    grid
    xlabel('Frequency'), ylabel('Power')
    title('\tau-Recurrence Rate Spectrum')
end
