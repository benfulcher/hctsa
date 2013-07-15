function out = ML_l1pwc_sweep_lambda(y,lambdar)
% ML sweep lambda steps
% Input signal y, and a lambda range lambdar (vector)
% Returns statistics on the fit across this range...
% Uses Max Little's step detection code
% Ben Fulcher 13/4/2010

Llambdar = length(lambdar);
nsegs = zeros(Llambdar,1);
rmserrs = zeros(Llambdar,1);
rmserrpsegs = zeros(Llambdar,1);

for i = 1:length(lambdar)
    lambda = lambdar(i);
    outi = ML_step_detection(y,'l1pwc',lambda);
    nsegs(i) = outi.nsegments;
    rmserrs(i) = outi.rmsoff;
    rmserrpsegs(i) = outi.rmsoffpstep;
end

% rmserrs gets under ** for first time
rmsunderx = @(x) find(rmserrs < x, 1, 'first');

out.rmserrsu05 = lambdar(rmsunderx(0.5));
if isempty(out.rmserrsu05), out.rmserru05 = NaN; end
out.rmserrsu02 = lambdar(rmsunderx(0.2));
if isempty(out.rmserrsu02), out.rmserru02 = NaN; end
out.rmserrsu01 = lambdar(rmsunderx(0.1));
if isempty(out.rmserrsu01), out.rmserru01 = NaN; end

% nsegs gets under ** for the first time
nsegunderx = @(x) find(nsegs < x, 1, 'first');

out.nsegsu005 = lambdar(nsegunderx(0.05));
if isempty(out.nsegsu005), out.nsegsu005 = NaN; end
out.nsegsu001 = lambdar(nsegunderx(0.01));
if isempty(out.nsegsu001), out.nsegsu001 = NaN; end

% correlation between #segments, rmserrs
R = corrcoef(nsegs,rmserrs);
out.corrsegerr = R(2,1);

% maximum rmserrpsegment
indbest = find(rmserrpsegs == max(rmserrpsegs),1,'first');
out.bestrmserrpseg = rmserrpsegs(indbest);
out.bestlambda = lambdar(indbest);


end