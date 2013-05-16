function test(mode)

% test nearest neighbor search based mex files

% recompile

error_flag = 0;
if nargin < 1
	mode = 'all';
end

disp('Fast nearest neighbor search routines test')


load points.dat
dat = points;
%dat = generate_chaotic_data(40000, 20);
%size(dat)

if strcmp(mode, 'delaunay2D') | strcmp(mode, 'all')
	error_flag = error_flag + against_delaunay(10000);
end
if strcmp(mode, 'brute') | strcmp(mode, 'all')
	error_flag = error_flag + against_brute(6000, 6, 1000);
end
if strcmp(mode, 'fullpartial') | strcmp(mode, 'all')
	error_flag = error_flag + against('fnearneigh', 'fnearneigh_partial', dat, 5, 5000, 200);
end
if strcmp(mode, 'fullpartial_brute') | strcmp(mode, 'all')
	error_flag = error_flag + against('nn_brute', 'nn_brute_partial', dat, 5, 400, 200);
end
if strcmp(mode, 'including') | strcmp(mode, 'all')
	error_flag = error_flag + including(dat, 32, 2000, 200, 'eucl');
end
if strcmp(mode, 'including_max') | strcmp(mode, 'all')
	error_flag = error_flag + including(dat, 16, 2000, 200, 'max');
end
if strcmp(mode, 'range') | strcmp(mode, 'all')
	error_flag = error_flag + range(dat, 1.5, 5000, 200, 'eucl');
end
if strcmp(mode, 'range_max') | strcmp(mode, 'all')
	error_flag = error_flag + range(dat, 1.5, 5000, 200, 'max');
end
if strcmp(mode, 'corr') | strcmp(mode, 'all')
	error_flag = error_flag + corr(dat, 1.5, 5000, 200, 'eucl');
end
if strcmp(mode, 'corr_max') | strcmp(mode, 'all')
	error_flag = error_flag + corr(dat, 1.5, 5000, 200, 'max');
end
disp(' ')
disp(' ')

if error_flag
	disp([num2str(error_flag) ' tests failed !!!'])
else
	disp(['All tests passed without errors'])
end

function error_flag = against(command1, command2, P, NNR, NREF, PAST)

error_flag = 0;
command_name1 = which(command1);
command_name2 = which(command2);
disp(['Testing ' command_name1 ' against ' command_name2])

ref = randref(1,size(P,1) ,NREF);

for DIM = 2:size(P,2)
    [i,d] = feval(command1, P(:,1:DIM), ref, NNR, PAST);
	[i2,d2] = feval(command2, P(:,1:DIM), ref, NNR, PAST);
	if all(all(i==i2))
		disp(['Returnded values are identical for dimension ' num2str(DIM)])
	else
		disp('ERROR : Returnded indices are not identical')
		error_flag = error_flag + 1;
	end
end

function error_flag = including(P, NNR, NREF, PAST, NORM)

error_flag = 0;

disp('Testing if repeated calls return same results')
disp(['Doing preprocessing for norm : ' NORM]);
disp(' ')

atria = nn_prepare(P, NORM);
ref = randref(1,size(P,1) ,NREF);

[nn, dist] = nn_search(P, atria, ref, 1, PAST);

for nnr=2:NNR
	[nn2, dist2] = nn_search_partial(P, atria, ref, nnr, PAST);
	
	if all(all(dist==dist2(:,1:end-1)))
		disp(['Returnded values are identical for NNR=' num2str(nnr)])
	else
		disp('ERROR : Returnded indices are not identical')
		error_flag = error_flag + 1;
	end	
	nn =  nn2;
	dist = dist2;
end

function error_flag = range(P, R, NREF, PAST, NORM)

error_flag = 0;

disp('Testing range search')
disp(['Doing preprocessing for norm : ' NORM]);
disp(' ')

atria = nn_prepare(P, NORM);
ref = randref(1,size(P,1) ,NREF);

tic; count = range_search(P, atria, ref, R, PAST); toc
tic; [count2, v] = range_search(P, atria, ref, R, PAST); toc
tic; [count3, v3] = range_search_partial(P, atria, ref, R, PAST); toc

mean(count)

if all(count==count2) & all(count==count3) 
	disp('Returnded values are identical, OK')
else
	disp('ERROR : Returnded indices are not identical')
	error_flag = 1;
end 

function error_flag = corr(P, R, NREF, PAST, NORM)

error_flag = 0;

disp('Testing corrsum with/out partial search')
disp(['Doing preprocessing for norm : ' NORM]);
disp(' ')

atria = nn_prepare(P, NORM);
ref = randref(1,size(P,1) ,NREF);

tic; [c,d, e] = corrsum(atria, P, ref, [R/1000 R], PAST); toc
tic; [c2, d2 ,e2] = corrsum_partial(atria, P, ref, [R/1000 R], PAST); toc

if all(c==c2) & all(e==e2) 
	disp('Returnded values are identical, OK')
else
	disp('ERROR : Returnded indices are not identical')
	error_flag = 1;
end 


function error_flag = against_delaunay(N)

error_flag = 0;

disp('Testing against matlab''s delaunay triangulation (only 2D, euclidian distance)')
disp(' ')

p = 20 * (rand(4 * N, 2)-0.5); 	% create 2-dimensional point set
ref = 25 * (rand(N, 2)-0.5); 	% create 2-dimensional reference point set (query points)

disp('Benchmarking Delaunay triangulation')

tic
tri = delaunay(p(:,1), p(:,2)); 
k = dsearch(p(:,1), p(:,2), tri, ref(:,1), ref(:,2)); 
toc

disp('Benchmarking fast nearest neighbor searcher (fnearneigh)')

tic
k2 = fnearneigh(p, ref, 1); 		
toc

if all(k==k2)
	disp('Returnded values are identical, OK')
else
	disp('ERROR : Returnded indices are not identical')
	error_flag = 1;
end

disp('Benchmarking fast nearest neighbor searcher with partial distance computation (fnearneigh_partial)')

tic
k2 = fnearneigh_partial(p, ref, 1); 		
toc

if all(k==k2)
	disp('Returnded values are identical, OK')
else
	disp('ERROR : Returnded indices are not identical')
	error_flag = error_flag + 1;
end

function error_flag = against_brute(N,D, Nref)

error_flag = 0;

disp('Testing against brute force search (eucl)')
disp(' ')

p = 20 * (rand(N, D)-0.5); 	% create D-dimensional point set
ref = randref(1,N, Nref);

disp('Benchmarking matlab brute force searcher')

tic
[nn, dist] = brute(p, ref, 5, 0);
toc

disp('Benchmarking mex brute force searcher (nn_brute)')

tic
[nn2,dist2] = nn_brute(p, ref, 5, 0);		
toc

if all(nn==nn2)
	disp('Returnded values are identical, OK')
else
	disp('ERROR : Returnded indices are not identical')
	error_flag = 1;
end

disp('Benchmarking mex partial distance brute force searcher (nn_brute_partial)')

tic
[nn2,dist2] = nn_brute_partial(p, ref, 5, 0);		
toc

if all(nn==nn2)
	disp('Returnded values are identical, OK')
else
	disp('ERROR : Returnded indices are not identical')
	error_flag = 1;
end



disp('Benchmarking fast nearest neighbor searcher (fnearneigh)')

tic
[nn2,dist2] = fnearneigh(p, ref, 5, 0);		
toc

if all(nn==nn2)
	disp('Returnded values are identical, OK')
else
	disp('ERROR : Returnded indices are not identical')
	error_flag = 1;
end


disp('Benchmarking fast nearest neighbor searcher with partial distance (fnearneigh_partial)')

tic
[nn2,dist2] = fnearneigh_partial(p, ref, 5, 0);		
toc

if all(nn==nn2)
	disp('Returnded values are identical, OK')
else
	disp('ERROR : Returnded indices are not identical')
	error_flag = 1;
end



function e = generate_chaotic_data(N, DIM)
 
x = chaosys(N + 1000, 0.03, [0.1 0.03298 0.001]);
e = data(embed(signal(x(1001:end,1)), DIM, 1));


