function [SSPE, REAL, PRED]  = ...
    mixed(ts_core,DIMS,future,LAGS,normhandle,SHIFTS,exclude,mode, ...
	   nn_number,exp)

%  MIXED calculates single step prediction error using mixed states.
%
%
%      MIXED(TS,DIMS) TS is a mxn matrix which holds the n
%      time series of m values. DIMS is a 1xn matrix which holds
%      the embedding dimensions of the time series. The return
%      value is the single step prediction error.
%      future is set to 1, LAGS is also set to 1 for each time
%      series, the norm function is the standard deviation (std()),
%      the SHIFTS are set to 0 and exclude is set to 0. The default
%      mode is normal and the default number of nn_number is 1.
%      The default for exp is 0.
%
%      MIXED(TS,DIMS,future) future is the number of time
%      steps that should be predicted ahead. It is an integer
%
%      MIXED(TS,DIMS,future,LAGS) LAGS is an 1xn matrix which
%      holds the time lags for each time series.
%
%      MIXED(TS,DIMS,future,LAGS,normhandle) normhandle is a
%      function handle that is used to normalize the data. The
%      function must be of the form f(row_vector) and must return a
%      scalar value. Useful values are max() or std(). A handle of
%      the value -1 indicates that no normalisation is desired.
%
%      MIXED(TS,DIMS,future,LAGS,normhandle,SHIFTS) SHIFTS is
%      a 1xn-1 matrix which defines the relative embedding shift of the
%      additional time series to the one that should be predicted.
%
%      MIXED(TS,DIMS,future,LAGS,normhandle,SHIFTS,exclude)
%      exclude is an integer value which gives the number of
%      neighbours in time space that should be excluded from the NN
%      search.
%
%      MIXED(TS,DIMS,future,LAGS,normhandle,SHIFTS,exclude,mode)
%      mode switches between some prediction flavours (see
%      paper). Possible values are normal, absolute and integrated.
%
%      MIXED(TS,DIMS,future,LAGS,normhandle,SHIFTS,exclude,mode,nn_number)
%      nn_number defines the number of nearest neighbours.
%
%      MIXED(TS,DIMS,future,LAGS,normhandle,SHIFTS,exclude,mode,nn_number,exp)
%      exp defines the exponent in the normalisation weights.
%
%      [sspe real pred] = MIXED(...) The return value sspe is
%      the single step prediction error. real and pred are matrices
%      of the real and predicted time series values respectively.
%
%  17.6.2002 Kevin Bube - Drittes Physikal. Institut Uni Goettingen
%
%  Copyright 1997-2002 DPI Goettingen
%  License http://www.physik3.gwdg.de/tstool/gpl.txt
%

%
% Changelog:
%             9.11.01 complete rewrite of mixed.m
%             17.6.02 checked into CVS
%

%
% BUGS:       poor help documentation
%             surely not every faulty user input is handled
%             default values only refer to two time series
%

%
% check preconditions and set some useful variables
%
narginchk(2,10)

TS = data (ts_core);

if ~isreal(TS) | ~isreal(DIMS)
  error ('matrices must be real')
end

switch nargin
 case 2
  future = 1;
  LAGS = [1 1];
  normhandle = @std;
  exclude = 0;
  SHIFTS = 0;
  mode = 'normal';
  nn_number = 1;
  exp = 0;
  SIZE_TS = size(TS);
  SIZE_DIMS = size(DIMS);
  SIZE_LAGS = size(LAGS);
  SIZE_SHIFTS = size(SHIFTS);

 case 3
  if max (size(future)) ~= 1 | future <= 0
    error ('future must be a positive integer')
  end
  LAGS = [1 1];
  normhandle = @std;
  SHIFTS = 0;
  exclude = 0;
  mode = 'normal';
  nn_number = 1;
  exp = 0;
  SIZE_TS = size(TS);
  SIZE_DIMS = size(DIMS);
  SIZE_LAGS = size(LAGS);
  SIZE_SHIFTS = size(SHIFTS);

 case 4
  SIZE_TS = size(TS);
  SIZE_DIMS = size(DIMS);
  SIZE_LAGS = size(LAGS);
  if max (size(future)) ~= 1 | future <= 0
    error ('future must be a positive integer')
  end
  normhandle = @std;
  SHIFTS = 0;
  exclude = 0;
  mode = 'normal';
  nn_number = 1;
  exp = 0;
  if (SIZE_LAGS(1) ~= 1 | SIZE_LAGS(2) ~= SIZE_TS(2))
    error ('LAGS dims are incorrect')
  end
  SIZE_SHIFTS = size(SHIFTS);

 case 5
  SIZE_TS = size(TS);
  SIZE_DIMS = size(DIMS);
  SIZE_LAGS = size(LAGS);
  if max (size(future)) ~= 1 | future <= 0
    error ('future must be a positive integer')
  end
  if max (size(normhandle)) ~= 1 | normhandle <= 0
    error ('normhandle must be a function pointer')
  end
  SHIFTS = 0;
  exclude = 0;
  mode = 'normal';
  nn_number = 1;
  exp = 0;
  if (SIZE_LAGS(1) ~= 1 | SIZE_LAGS(2) ~= SIZE_TS(2))
    error ('LAGS dims are incorrect')
  end
  SIZE_SHIFTS = size(SHIFTS);

 case 6
  SIZE_TS = size(TS);
  SIZE_DIMS = size(DIMS);
  SIZE_LAGS = size(LAGS);
  if max (size(future)) ~= 1 | future <= 0
    error ('future must be a positive integer')
  end
  if max (size(normhandle)) ~= 1
    error ('normhandle must be a function pointer')
  end
  SIZE_SHIFTS = size(SHIFTS);
  if (SIZE_SHIFTS(1) ~= 1 | SIZE_SHIFTS(2) ~= SIZE_TS(2)-1)
    error ('SHIFT dims are incorrect')
  end
  exclude = 0;
  mode = 'normal';
  nn_number = 1;
  exp = 0;
  if (SIZE_LAGS(1) ~= 1 | SIZE_LAGS(2) ~= SIZE_TS(2))
    error ('LAGS dims are incorrect')
  end

 case 7
  SIZE_TS = size(TS);
  SIZE_DIMS = size(DIMS);
  SIZE_LAGS = size(LAGS);
  if max (size(exclude)) ~= 1 | exclude < 0
    error ('exclude must be a nonegative integer')
  end
  if max (size(future)) ~= 1 | future <= 0
    error ('future must be a positive integer')
  end
  if max (size(normhandle)) ~= 1
    error ('normhandle must be a function pointer')
  end
  SIZE_SHIFTS = size(SHIFTS);
  if (SIZE_SHIFTS(1) ~= 1 | SIZE_SHIFTS(2) ~= SIZE_TS(2)-1)
    error ('SHIFT dims are incorrect')
  end
  mode = 'normal';
  nn_number = 1;
  exp = 0;
  if (SIZE_LAGS(1) ~= 1 | SIZE_LAGS(2) ~= SIZE_TS(2))
    error ('LAGS dims are incorrect')
  end

 case 8
  SIZE_TS = size(TS);
  SIZE_DIMS = size(DIMS);
  SIZE_LAGS = size(LAGS);
  if max (size(exclude)) ~= 1 | exclude < 0
    error ('exclude must be a nonegative integer')
  end
  if max (size(future)) ~= 1 | future <= 0
    error ('future must be a positive integer')
  end
  if max (size(normhandle)) ~= 1
    error ('normhandle must be a function pointer')
  end
  SIZE_SHIFTS = size(SHIFTS);
  if (SIZE_SHIFTS(1) ~= 1 | SIZE_SHIFTS(2) ~= SIZE_TS(2)-1)
    error ('SHIFT dims are incorrect')
  end
  if (SIZE_LAGS(1) ~= 1 | SIZE_LAGS(2) ~= SIZE_TS(2))
    error ('LAGS dims are incorrect')
  end
  nn_number = 1;
  exp = 0;

 case 9
  SIZE_TS = size(TS);
  SIZE_DIMS = size(DIMS);
  SIZE_LAGS = size(LAGS);
  if max (size(exclude)) ~= 1 | exclude < 0
    error ('exclude must be a nonegative integer')
  end
  if max (size(future)) ~= 1 | future <= 0
    error ('future must be a positive integer')
  end
  if max (size(normhandle)) ~= 1
    error ('normhandle must be a function pointer')
  end
  SIZE_SHIFTS = size(SHIFTS);
  if (SIZE_SHIFTS(1) ~= 1 | SIZE_SHIFTS(2) ~= SIZE_TS(2)-1)
    error ('SHIFT dims are incorrect')
  end
  if (SIZE_LAGS(1) ~= 1 | SIZE_LAGS(2) ~= SIZE_TS(2))
    error ('LAGS dims are incorrect')
  end
  if max (size(nn_number)) ~= 1 | nn_number <= 0
    error ('nn_number must be a positive integer')
  end
  exp = 0;

 case 10
  SIZE_TS = size(TS);
  SIZE_DIMS = size(DIMS);
  SIZE_LAGS = size(LAGS);
  if max (size(exclude)) ~= 1 | exclude < 0
    error ('exclude must be a nonegative integer')
  end
  if max (size(future)) ~= 1 | future <= 0
    error ('future must be a positive integer')
  end
  if max (size(normhandle)) ~= 1
    error ('normhandle must be a function pointer')
  end
  SIZE_SHIFTS = size(SHIFTS);
  if (SIZE_SHIFTS(1) ~= 1 | SIZE_SHIFTS(2) ~= SIZE_TS(2)-1)
    error ('SHIFT dims are incorrect')
  end
  if (SIZE_LAGS(1) ~= 1 | SIZE_LAGS(2) ~= SIZE_TS(2))
    error ('LAGS dims are incorrect')
  end
  if max (size(nn_number)) ~= 1 | nn_number <= 0
    error ('nn_number must be a positive integer')
  end
  if max (size(exp)) ~= 1 | exp <= 0
    error ('exp must be a positive integer')
  end

end


if (SIZE_DIMS(1) ~= 1 | SIZE_DIMS(2) ~= SIZE_TS(2))
  error ('DIMS dimensions are incorrect')
end


%
% check if we have to normalize the data
%
if ~isa(normhandle,'function_handle')     % -1 is an invalid function
                                          % handle which means we don't
                                          % normalize
					  for k = 1:SIZE_TS(2)
					    TS(:,k) = TS(:,k) ./ feval(normhandle, TS(:,k));
					  end
end

%
% First we treat the special case that all embedding dimensions are
% equal to zero, meaning the mixed state vector would be empty.
%
if (max (DIMS) == 0)

  % create random indizes
  indices = ceil (rand(SIZE_TS(1),1) * (SIZE_TS(1)-future));

  % calculate sspe with randomly chosen points
  SSPE = mean (abs (TS((1 : SIZE_TS(1)-future) + future, 1) - ...
		    TS(indices(1 : SIZE_TS-future) + future, 1)));
  REAL = TS((1 : SIZE_TS(1)-future) + future, 1);
  PRED = TS(indices(1 : SIZE_TS-future) + future, 1);
  return
end


%
% That was quite simple. But now we treat the normal case(s)...
% We have to begin with the construction of our embedding vector p.
%
[p startIndex] = mixembed (TS, DIMS, LAGS, SHIFTS);
sP = size (p);

%
% As we have done the embedding, we can do the most important thing:
% the next neighbour search.
%
atria = nn_prepare (p(1:end-future, :));

%
% Here we distiguish between some prediction flavours an we
% determine the real and the predicted values of the time
% series. In order to do so we need some index magic...
%

if ~strcmp (mode, 'normal')
  % is needed because we need nn_number+1 neighbours for the other
  % modes
  nn_number = nn_number + 1;
end

[index dist] = nn_search (p(1:end-future, :), atria, 1:sP(1)-future, ...
			  nn_number, exclude);
REAL = TS((startIndex:startIndex-1+sP(1)-future) + future, 1);


switch mode
 case 'normal'
  PRED = TS(index(1:end) + startIndex-1 + future, 1);

 case 'absolute'
  size_dist = size (dist);
  zeros (size_dist(1),size_dist(2)-1);
  for k = 1:size_dist(1)
    r(k,:) = dist(k, 1:end-1) / dist(k, end);
  end

  weights      = (1 - r.^exp).^exp;
  norm         = sum (weights');
  norm         = norm';
  size_weights = size (weights);

  sumbuffer = zeros(size(weights));
  for k = 1:size_weights(1)
    for l = 1:size_weights(2)
      sumbuffer(k,l) = weights(k,l) * TS(index(k,l) + startIndex-1 ...
					 + future, 1);

    end
  end

  sumbuffer  = sum (sumbuffer');
  PRED       = sumbuffer' ./ norm;


 case 'integrated'
  size_dist = size (dist);
  zeros (size_dist(1),size_dist(2)-1);
  for k = 1:size_dist(1)
    r(k,:) = dist(k, 1:end-1) / dist(k, end);
  end

  weights      = (1 - r.^exp).^exp;
  norm         = sum (weights');
  norm         = norm';
  size_weights = size (weights);

  sumbuffer = zeros(size(weights));
  for k = 1:size_weights(1)
    for l = 1:size_weights(2)
      sumbuffer(k,l) = (weights(k,l) ...
			* (TS(index(k,l) + startIndex-1 + future, 1) ...
			   - TS(index(k,l) + startIndex-1, 1)));
    end
  end

  sumbuffer  = sum (sumbuffer');
  PRED       = ((TS(startIndex:size_weights(1)+startIndex-1))' ...
		+ (sumbuffer' ./ norm));


 otherwise
  error ('unknown mode');
end

%
% Yuhuu, gotcha!
%
SSPE = mean (abs (REAL - PRED));
return
