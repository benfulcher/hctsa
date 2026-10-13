function tests = Robust_salv_enrandTests
% Tests for EN_RandomizeOrdinal (exact expected ordinal-pattern decay under randomization)
tests = functiontests(localfunctions);
end

function testDeterministicAndOrdinal(testCase)
	% no randomness, and invariant to a strictly increasing transformation
	t = (1:300)';
	y = sin(0.2 * t) + 0.3 * cos(0.05 * t.^1.1);
	o1 = EN_RandomizeOrdinal(y);
	o2 = EN_RandomizeOrdinal(y);
	o3 = EN_RandomizeOrdinal(exp(2 * y) + 5);
	verifyEqual(testCase, o1.halftime, o2.halftime);
	verifyEqual(testCase, o1.halftime, o3.halftime, 'AbsTol', 1e-12);
	verifyEqual(testCase, o1.normPermEnInf, 1, 'AbsTol', 1e-3); % near 1 (draws can coincide)
end

function testReferenceValueWithTies(testCase)
	% reference value from an independent implementation (repeated values included)
	t = (1:200)';
	y = round(3 * sin(0.3 * t) + mod(t * 7, 5) * 0.5);
	o = EN_RandomizeOrdinal(y);
	verifyEqual(testCase, o.halftime, 0.120684956226594, 'AbsTol', 1e-9);
end

function testBruteForceDistributions(testCase)
	% the one-replaced and fully randomized pattern distributions, by enumerating
	% every draw on a short series with repeated values
	y = [2; 1; 2; 3; 1; 2; 2];
	N = length(y);
	x = BF_Embed(y, 1, 3, false);
	Prep = zeros(1, 6); Pinf = zeros(1, 6);
	for i = 1:size(x, 1)
		for c = 1:3
			for j = 1:N
				z = x(i, :); z(c) = y(j);
				p = BF_OrdinalPatternRank(z);
				Prep(p) = Prep(p) + 1 / (3 * N);
			end
		end
	end
	Prep = Prep / size(x, 1);
	for a = 1:N, for b = 1:N, for c = 1:N
		p = BF_OrdinalPatternRank(y([a, b, c])');
		Pinf(p) = Pinf(p) + 1 / N^3;
	end, end, end
	H = @(P) -sum(P(P > 0) .* log2(P(P > 0))) / log2(6);
	o = EN_RandomizeOrdinal(y);
	verifyEqual(testCase, o.normPermEnRep1, H(Prep), 'AbsTol', 1e-12);
	verifyEqual(testCase, o.normPermEnInf, H(Pinf), 'AbsTol', 1e-12);
	pe = EN_PermEn(y, 3, 1);
	verifyEqual(testCase, o.normPermEn0, pe.normPermEn, 'AbsTol', 1e-12);
end

function testNoStructure(testCase)
	% a constant series has no ordinal structure to lose
	o = EN_RandomizeOrdinal(ones(50, 1));
	verifyTrue(testCase, isnan(o.halftime));
end
