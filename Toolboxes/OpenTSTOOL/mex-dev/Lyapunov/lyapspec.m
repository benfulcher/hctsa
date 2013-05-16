function lspec = lyapspec(ts, dim, delay, NNR, eps, Nref)

ts = ts(:);
N = length(ts) - dim*delay;
ref = 1:delay:min(N,(1+(Nref-1)*delay));

[nn, dists, points] = emb_nn_search(ts(1:(end-delay)), [dim delay], ref, NNR+1, -1);

Q = diag(ones(dim,1));	
lspec = [];zeros(dim, 1);

for i=1:length(ref)
	X = points(nn(i,:), :);
	X0 = mean(X);		% center of gravity	
	Y = X - repmat(X0, NNR+1, 1);
	[U,W,V] = svd(Y);
	w = diag(W) / sum(diag(W));
	k = max(find(w >= eps));
	V = V(:,1:k);
	Z = Y * V;
	
	images = ts(dim*delay+nn(i,:));
	%all(points(nn(i,:)+delay,1) == images);
	
	
	a = [Z ones(NNR+1,1)] \ images;
	
	b = a(end);
	a = a(1:end-1);
	c = V * a;

	%mean(abs(images-(X * c + ( b - X0 * c))) ./ abs(images))
	
	JAC = sparse([ones(1,dim) 2:dim], [1:dim 1:dim-1], [c' ones(1, dim-1)]);
	%full(JAC);
	
	A =  JAC * Q;
	[Q,R,E] = qr(A);

	lspec = [lspec; diag(R)'];
end

lspec = mean(log(abs(lspec(min(length(ref)/10,500):end,:))));
