function g = gradf(x)
n = length(x);
f = sqrt(n)^n*prod(x);
g = zeros(n,1);
for ni = 1:n
    Prodx = x; Prodx(ni) = 1;
    g(ni) = sqrt(n)^n*prod(Prodx);
end
g = -g;
