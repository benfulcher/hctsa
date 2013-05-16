function r = mrdivide(p,q)

p = unit(p);
q = unit(q);

r = unit([p.factor / q.factor (p.exponents - q.exponents)]);
