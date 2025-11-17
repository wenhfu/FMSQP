function c = func(x)
n = length(x);
c1 = x(1)^2-x(2)+1;
c2 = 1-x(1)+(x(2)-4)^2;
lb = [0;0];
ub = [10;10];
c = [c1;c2;-x+lb;x-ub];
