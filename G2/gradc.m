function gc = gradc(x)
n = length(x);
m = 1 + 1 + 2*n;

gc = zeros(n,m);
gc1 = zeros(n,1);
for ni = 1:n
    Prodx = x; Prodx(ni) = 1;
    gc1(ni) = -prod(Prodx);
end
% gc1 = -prod(x)./x;
gc2 = ones(n,1);
gc3 = -eye(n);
gc4 = eye(n);
gc = [gc1,gc2,gc3,gc4];
