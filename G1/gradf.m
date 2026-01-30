function g = gradf(x)
g1 = [5*ones(4,1);zeros(9,1)];
g2 = [-5*2*x(1:4);zeros(9,1)];
g3 = [zeros(4,1);-ones(9,1)];
g = g1 + g2 + g3;
