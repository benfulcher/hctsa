function varargout=recons(varargin)
%RECONS   Reconstruct a time series from a recurrence plot.
%    Y = RECONS(X) reconstructs a time series Y from the 
%    recurrence plot in the matrix X.
%
%    Y = RECONS(X,NAME) reconstructs the time series using
%    the named cumulative distribution function, which can
%    be 'norm' or 'Normal' (defeault), 'unif' or 'Uniform'.
%
%    Y = RECONS(X,P) reconstructs the time series using
%    the cumulative distribution function given by vector P.
%
%    See also CRP, CRP2, JRP, TWINSURR.
%
%    References: 
%    Thiel, M., Romano, M. C., Kurths, J.: 
%    How much information is contained in a recurrence plot?, 
%    Phys. Lett. A, 330, 2004.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2005-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:34:10 $
% $Revision: 5.3 $
%
% $Log: recons.m,v $
% Revision 5.3  2009/03/24 08:34:10  marwan
% copyright address changed
%
% Revision 5.2  2008/06/24 14:11:55  marwan
% Including of a abort condition in order to avoid infinite loops.
%
% Revision 5.1  2008/01/25 12:47:25  marwan
% initial import
%
%
%
%
% This program is part of the new generation XXII series.
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


errcode = 0; %xout = [];
error(nargchk(1,2,nargin));
if nargout>1, error('Too many output arguments'), end

if isnumeric(varargin{1})==1 		% read RP
    X = double(varargin{1});
else 
    error('Input must be numeric.')
end



N = size(X); maxN = prod(N)+1;
if N(1) ~= N(2)
    error('Input must be square matrix.')
end

try
% final distribution        
errcode = 1;
i=(-5:10/(N(1)-1):5)';
fs=cumsum(10*(1/(sqrt(2*pi))).*exp(-(i).^2./2));
P=interp1(fs,i,[1:length(fs)],'nearest')';

if nargin > 1 & isnumeric(varargin{2})==1
    P = varargin{2};
elseif nargin > 1 & ischar(varargin{2})==1
    if strcmpi(varargin{2}(1:3), 'uni')
        P = (0:1/(N(1)-1):1)';
        disp('uni')
    end
end
P(isnan(P)) = []; % remove NaN from distribution


h_waitbar = waitbar(0,'Remove double columns'); drawnow

% remove double rows
errcode = 2;
[X2 uniI J]=unique(X,'rows');
uniI = sort(uniI);

% this are the removed rows
rem = setdiff(1:length(X2),uniI); 
%X3(rem,:) = X(rem,:);

% this are the corresponding columns
for i = 1:length(rem), waitbar(i/length(rem))
%    [dummy, dblI(i), a2] = intersect(X, X(rem(i),:), 'rows');
    X2 = X;
    X2(rem,:) = -1;
    dblI(i) = find(all((X2 == repmat(X(rem(i),:), size(X2,2), 1))'));
end

% number of non-neighbours
errcode = 3;

waitbar(0, h_waitbar, 'Search neighbours'); drawnow
NN = maxN*ones(N);
%N = zeros(size(X1));
for i = 1:size(X,1), waitbar(i/size(X,1))
    if ~ismember(i,rem)
        % indices of RPs at i
        I = find(X(i,:));
        I(ismember(I,rem)) = [];

        % number of non-neighbours
        n = sum((X(I,:) - repmat(X(i,:),length(I),1)) == 1,2);
        NN(i,I) = n';
    end
    
end
errcode = 4;

% remove entries on the main diagonal in N
i(1:size(NN,1)) = 1; NN(find(diag(i))) = maxN;

% determine the both indices where N == 0 for all i
iNaN = NN == maxN; i = find(iNaN); NN2 = NN; NN2(i) = zeros(size(i)); % remove NaNs
n = sum(NN2,1);
borderI=find(n==0);
borderI(ismember(borderI,rem)) = [];

k = borderI(1);
r = k; % start point in the rank order vector

errcode = 5;
waitbar(0, h_waitbar, 'Reconstructing time series'); drawnow
while min(NN(:)) < maxN
    % look for the minimal N(k,:)
    waitbar(length(r)/N(1))
    [dummy kneu] = min(NN(k,:));
    if ~isempty(dummy) kneu = find(NN(k,:) == dummy); else break, end
    if length(kneu) > 1 % if a set of minimal N exist
        [dummy kneu_index] = min(NN(kneu,k));
        kneu=kneu(kneu_index);
    end
    NN(:,k) = maxN;
    kold = k;
    k=kneu;
    r = [r; k]; % add the new found index to the rank order vector
    if length(r) > N(1) + 1
        delete(h_waitbar)
        disp(['Critical abort. Could not find enough corresponding neigbours.',10,'Perhaps the recurrence threshold is too small.'])
        if nargout
            varargout{1} = NaN;
        end
        return
    end
    
end

errcode = 6;
% adjust length of the distribution
P = interp1((1:length(P))/length(P), P, (1:length(r))/length(r),'cubic');

% assign the distribution to the rank order series
waitbar(0, h_waitbar, 'Assigning distribution'); drawnow
clear xneu
for i = 1:length(r);
    xneu(r(i)) = P(i);
    waitbar(i/length(r))
end

errcode = 7;

% add the removed (doubled) elements
for i =1:length(rem)
    xneu(rem(i)) = xneu(dblI(i));
end
delete(h_waitbar)


errcode = 8;

if nargout
    varargout{1} = xneu;
else
    xneu
end




%%%%%%% error handling

%if 0
catch
  z=whos;x=lasterr;y=lastwarn;in=varargin{1};
  if ~isempty(findobj('Tag','TMWWaitbar')), delete(findobj('Tag','TMWWaitbar')), end
  if ~strcmpi(lasterr,'Interrupt')
    fid=fopen('error.log','w');
    err=fprintf(fid,'%s\n','Please send us the following error report. Provide a brief');
    err=fprintf(fid,'%s\n','description of what you were doing when this problem occurred.');
    err=fprintf(fid,'%s\n','E-mail or FAX this information to us at:');
    err=fprintf(fid,'%s\n','    E-mail:  marwan@pik-potsdam.de');
    err=fprintf(fid,'%s\n','       Fax:  ++49 +331 288 2640');
    err=fprintf(fid,'%s\n\n\n','Thank you for your assistance.');
    err=fprintf(fid,'%s\n',repmat('-',50,1));
    err=fprintf(fid,'%s\n',datestr(now,0));
    err=fprintf(fid,'%s\n',['Matlab ',char(version),' on ',computer]);
    err=fprintf(fid,'%s\n',repmat('-',50,1));
    err=fprintf(fid,'%s\n',x);
    err=fprintf(fid,'%s\n',y);
    err=fprintf(fid,'%s\n',[' during ==> recons']);
    if ~isempty(errcode)
      err=fprintf(fid,'%s\n',[' errorcode ==> ',num2str(errcode)]);
    else
      err=fprintf(fid,'%s\n',[' errorcode ==> no errorcode available']);
    end
    err=fclose(fid);
    disp('----------------------------');
    disp('       ERROR OCCURED ');
    disp(['   during executing recons']);
    disp('----------------------------');
    disp(x);
    if ~isempty(errcode)
      disp(['   errorcode is ',num2str(errcode)]);
    else
      disp('   no errorcode available');
    end
    disp('----------------------------');
    disp('   Please send us the error report. For your convenience, ')
    disp('   this information has been recorded in: ')
    disp(['   ',fullfile(pwd,'error.log')]), disp(' ')
    disp('   Provide a brief description of what you were doing when ')
    disp('   this problem occurred.'), disp(' ')
    disp('   E-mail or FAX this information to us at:')
    disp('       E-mail:  marwan@pik-potsdam.de')
    disp('          Fax:  ++49 +331 288 2640'), disp(' ')
    disp('   Thank you for your assistance.')
    warning('on')
  end
  try, set(0,props.root), delete(h_waitbar), end
  set(0,'ShowHidden','Off')
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
