function g = gradf(x)
n = length(x);
g = zeros(n,1);
g(1) = 37.293239 + 0.8356891*x(5);
g(3) = 2*5.3578547*x(3);
g(5) = 0.8356891*x(1);
