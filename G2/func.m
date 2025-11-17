function c = func(x)
n = length(x);
c1 = 0.75 - prod(x);
c2 = sum(x) - 7.5*n;
c3 = -x;
c4 = x-10;
c = [c1;c2;c3;c4];