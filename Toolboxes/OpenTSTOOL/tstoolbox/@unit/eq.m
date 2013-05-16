function r = eq(a,b)

if (a.factor == b.factor) & (a.exponents == b.exponents)
	r = 1;
else
	r = 0;
end
