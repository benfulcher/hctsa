function [Z,A,B,rss] = arch2(X, na, mode, sil)

%   [Z,A,B,rss] = arch2(X, na, mode, silent )
%
%   archetypal analysis of column orientated data set <X>
%
%   input arguments :
%
%   - each column of data is one 'observation', e.g. the sample values of
%     all channels in a multichannel measurement at one point in time
%
%   - na : number of generated archetypes
%
%   - mode can be one of the following : 'normalized' (default),
%     'mean', 'raw'
%     - in mode 'normalized' each column of data is centered by removing its mean
%       and then normalized by dividing through its standard deviation before
%       the covariance matrix is calculated
%     - in mode 'mean' only the mean of every column of data is removed
%     - in mode 'raw' no preprocessing is applied to data
%     - in mode 'scale' <X> is divided by max(abs(X))
%
%   - silent is an optional flag which supresses output of text and plot on the matlab
%     screen. Returned values (see below) are in no way affected
%
%
%   output arguments :
%
%   - Z : each column of Z is an archetype
%
%   - A : the columns of A are the coefficients of the archetypes to
%         the constrains ||X-Z*A|| -> min
%
%   - B : the columns of B create the X-mixtures (the archtypes)
%         Z=X*B
%
%   - rss : the residual sum of squares for each iteration
%
%   Christian Merkwirth & Joerg Wichard
%   Februar 1998


global silent

narginchk(2,3);

if nargin < 3
    mode = 'normalized';
end

if nargin < 4
    silent = 0;
else
    silent = 1;
end

[m,n] = size(X);

printline('archetypal analysis')
printline(['on data set of size ' num2str(n) 'x' num2str(m)]);

if (na < 1)
  printline('number of archetypes must be greater than zero');
end

if strncmp(mode, 'r',1)
  mode = 'raw';
  printline('no data preprocessing');
elseif strncmp(mode, 'm',1)
  mode = 'mean';
  printline('removing mean from data set');
  mn = mean(X,2);
  X = X - repmat(mn, 1, n);
elseif strncmp(mode, 's',1)
  mode = 'scale';
  printline('scaling data set');
  xm = max(max(abs(X)));
  X = (1/xm)*X;
else
  mode = 'normalized';
  printline('removing mean and normalizing data');
  mn = mean(X,2);
  X = X - repmat(mn, 1, n);
  dv = std(X,0,2);
  X = X ./ repmat(dv, 1, n);
end

gew1=20*m;              % Gewichtung der Convexit�tsbedingung
gew2=5*m;               % Gewichtung der Convexit�tsbedingung
tol=0.01;              %%Toleranz f�r die Abbruchbedingung
numb = 20;              %% Max. Anzahl der Iterationen

%%  x zuf�llig ausw�len
B=eye(n);
rp=randperm(na);
B=B(:,rp);
Z=X*B;

%% Hier beginnt die Alternierende Optimierung
rss(1,:)=[ 0 sum(sum(X .* X))]
count=0;

for c=1:numb;

  for l=1:na;

    %% Maximum an Z anf�gen
    MX=max(max(Z));
    Z(m+1,:)=gew1*MX;
    X(m+1,:)=gew1*MX;

    %% A suchen bei konstantem Z
    for i=1:n;
      A(:,i)= lsqnonneg(Z,X(:,i));
    end;

    %% Gewichtung entfernen
    X=X(1:m,:);
    Z=Z(1:m,:);

    %% Z suchen bei konstantem A
    for i=1:n;
      V(:,i)=A(l,i)*(X(:,i) - Z*A(:,i) + A(l,i)*Z(:,l));
    end;

    %% Singul�re Archetypen durch max ||X-Z*A|| ersetzen
    if ( (sum(A(l,:)) .* sum(A(l,:))) ==0)
      VT=sum((X-Z*A).*(X-Z*A));
      [VTC,VTI]=max(VT);
      B(:,l)=0;
      B(VTI,l)=1;

    else
      VS=(sum(A(l,:).*A(l,:)))*sum(V,2);
      MV=max(max(X));
      X(m+1,:)=gew2*MV;
      VS(m+1)=gew2*MV;

      B(:,l)=nnls(X,VS, 4 * max(size(X)) * norm(X,1) * eps);

      %% Gewichtung entfernen
      X=X(1:m,:);
      VS=VS(1:m);
    end;
  end;

  Z=X*B;

  %% norm RSS berechnen
  rss((c+1),:)=[ c sum(sum((X-Z*A).*(X-Z*A))) ];

  if (abs(rss(c+1,2)-rss(c,2)) < tol*rss(c,2) )
    break;
  end;

  count=count+1;
end;


%% Ergebnis auf Konsistenz �berpr�fen
[C,I]=max(rss);
if( (I < count) & (rss(I) < rss(count)) )
  printline('Minimum not reached')
end

if( count == numb )
  printline('Maximum Number of Iterations')
end

if silent
  figure(1)
  subplot(4,1,1)
  plot(rss(:,1), rss(:,2))
  title('Residual Sum of Squares')

  subplot(4,1,2)
  plot(X)
  title('Input Data');

  subplot(4,1,3)
  plot((Z),'k')
  title('Archetypes');

  subplot(4,1,4)
  plot((Z),'k')
  hold on;
  plot(X)
  hold off;
  title('Archetypes & Input Data');
end


function printline(string)
global silent
if silent~=1
        disp(string)
end
