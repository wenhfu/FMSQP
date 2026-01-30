function g = gradf(x)
n = length(x);
g = zeros(n,1);
for ni = 1:n
    Prodx = x; Prodx(ni) = 1;
    g(ni) = prod(Prodx)*exp(prod(x));
end
