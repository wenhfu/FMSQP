function gc = gradc(x)
x1 = x(1); x2 = x(2);
gc1 = pi*exp(cos(2*pi*x1)/2 + cos(2*pi*x2)/2)*sin(2*pi*x1) + (2*x1*exp(-(x1^2/2 + x2^2/2)^(1/2)/5))/(x1^2/2 + x2^2/2)^(1/2);
gc2 = pi*exp(cos(2*pi*x1)/2 + cos(2*pi*x2)/2)*sin(2*pi*x2) + (2*x2*exp(-(x1^2/2 + x2^2/2)^(1/2)/5))/(x1^2/2 + x2^2/2)^(1/2);
gc = [gc1;gc2];
