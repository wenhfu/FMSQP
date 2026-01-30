function gc = gradc(x)
n = length(x);
gc = [-eye(n),eye(n)];
gh = gradh(x);
gc = [gc,gh,-gh];