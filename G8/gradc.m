function gc = gradc(x)
n = length(x);
gc1 = [2*x(1),-1];
gc2 = [-1,2*(x(2)-4)];
Dc = [gc1;gc2;-eye(n);eye(n)];
gc = Dc';