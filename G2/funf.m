function f = funf(x)
n = length(x);
eqn1 = sum(cos(x).^4) - 2*prod(cos(x).^2);
eqn2 = sum((1:n)'.*(x.^2));
f = -eqn1^2/eqn2;


% f = -eqn1/sqrt(eqn2);