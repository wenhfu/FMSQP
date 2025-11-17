function h = funh(x)
h1 = sum(x.^2)-10;
h2 = x(2)*x(3)-5*x(4)*x(5);
h3 = x(1)^3+x(2)^3+1;
h = [h1;h2;h3];
