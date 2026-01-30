function c = func(x)
c1 = x(3) - x(4) - 0.55;
c2 = -x(3) + x(4) - 0.55;
lb = [0;0;-0.55;-0.55];
ub = [1200;1200;0.55;0.55];
c = [c1;c2;-x+lb;x-ub];
h = funh(x);
c = [c;h;-h];
