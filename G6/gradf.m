function g = gradf(x)
n = length(x);
g = zeros(n,1);
g(1) = 3*(x(1) - 10)^2;
g(2) = 3*(x(2) - 20)^2;
