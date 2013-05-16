function varargout = rtspec(varargin)
% RTSPEC   Recurrence time spectrum.
%    RTSPEC(X,M,T,E,FS,...) calculates the recurrence time spectrum
%    based on a recurrence plot using embedding dimension M,
%    embedding delay T, recurrence threshold E, and sampling
%    frequency FS. The input arguments are similar to those of the 
%    command CRP.
%
%    P = RTSPEC(...) returns the recurrence time spectrum
%    in vector P.
%
%    [P F] = RTSPEC(...) returns the recurrence time spectrum
%    in vector P and the vector of corresponding frequencies F.
%
%    Example: fs = 22;
%             x = sin(2*pi * [0:1/fs:44]);
%             rtspec(x,2,1,.1,fs)
%
%    See also CRP, RRSPEC.
%    

% Copyright (c) 2008 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2008/07/02 12:01:13 $
% $Revision: 5.1 $
%
% $Log: rtspec.m,v $
% Revision 5.1  2008/07/02 12:01:13  marwan
% initial import
%
% Revision 5.1  2008/07/01 13:09:27  marwan
% initial import
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% programme properties

global errcode props

init_properties
fs_init = 1;
w_init = 100;
sil = 1;

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
         method = 'max'; norm_str = 'nor';
    end

    % get the parameters for creating RP
    if max(size(varargin{1}))<=3
        disp('Error using ==> rtspec')
        disp('To less values in data X.')
        return
    end
    x=double(varargin{1});

    if ~isempty(varargin{i_double(2)}), m=varargin{i_double(2)}(1); else m=1; end
    if ~isempty(varargin{i_double(3)}), t=varargin{i_double(3)}(1); else t=1; end
    if ~isempty(varargin{i_double(4)}), e=varargin{i_double(4)}(1); else e=.1; end
    if ~isempty(varargin{i_double(5)}), w=varargin{i_double(5)}(1); else w=w_init; end
    if ~isempty(varargin{i_double(5)}), fs=varargin{i_double(5)}(1); else fs=fs_init; end
else
    disp('Error using ==> rtspec')
    disp('No valid arguments.')
    return
end


if size(x,1) < size(x,2), x = x'; end
N = length(x)-(m-1)*t;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% make RT spectrum

%% RP
nogui_str = 'nog';
if sil, nogui_str = 'sil'; end

X = crp2(x,m,t,e,method,norm_str,nogui_str);

%% recurrence times
[dummy1 dummy2 RT] = tt(X);

%% spectrum
    
f1 = 0:1/(1000*fs):1;
P1 = histc(1./(RT+1),f1);

P = histc(RT,1:2*max(RT));
f = [1:2*max(RT)]+.5;
f2 = fs./(f+1);

if nargout == 1
    varargout{1} = P(:);
elseif nargout == 2
    varargout{1} = P(:);
    varargout{2} = f2(:);
else
    
    semilogy(f1*fs,P1+1,fs./(f+1),P+1)
%    semilogy(fs./f,P+1)
    grid
    ylabel('Power')
    xlabel('Frequency')
    title('Recurrence Time Spectrum')

end
