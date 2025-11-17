function f = funf(x)
n = length(x);
f = sqrt(n)^n*prod(x);
f = -f;
