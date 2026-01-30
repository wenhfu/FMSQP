function c = func(x)
lb = [-1;-1];
ub = [1;1];
c = [-x+lb;x-ub];
h = funh(x);
c = [c;h;-h];
