function g = gradf(x)
n = length(x);
eqn1 = sum(cos(x).^4) - 2*prod(cos(x).^2);
eqn2 = sum((1:n)'.*(x.^2));
% f1 = eqn1^2;
% f2 = eqn2;
% 
% g = zeros(n,1);
% for ni = 1:n
%     df1 = 2*eqn1*(-4*cos(x(ni))^3*sin(x(ni)) + 4*prod(cos(x).^2)*tan(x(ni)));
%     df2 = 2*ni.*x(ni);
%     g(ni) = (df1*f2 - df2*f1)/(f2^2);
% end
% g = -g;
% 
% 
f1 = eqn1;
f2 = sqrt(eqn2);

g = zeros(n,1);
for ni = 1:n
    df1 = -4*cos(x(ni))^3*sin(x(ni)) + 4*prod(cos(x).^2)*tan(x(ni));
    df2 = ni*x(ni)/sqrt(eqn2);
    g(ni) = (df1*f2 - df2*f1)/(f2^2);
end
g = -g;


