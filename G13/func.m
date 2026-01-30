function c = func(x)
lb = [-2.3;-2.3;-3.2;-3.2;-3.2];
ub = [2.3;2.3;3.2;3.2;3.2];
c = [-x+lb;x-ub];
h = funh(x);
c = [c;h;-h];
