function f=funf(x)
f1 = sin(2*pi*x(1))^3*sin(2*pi*x(2));
f2 = x(1)^3*(x(1)+x(2));
f = f1/f2;
f = -f;

