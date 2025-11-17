function gc = gradc(x)
n = length(x);
gc1 = [0,0,1,-1];
gc2 = [0,0,-1,1];
Dc = [gc1;gc2;-eye(n);eye(n)];
gc = Dc';
gh = gradh(x);
gc = [gc,gh,-gh];