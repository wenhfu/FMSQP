function gh = gradh(x)
gh1 = 2*x';
gh2 = [0,x(3),x(2),-5*x(5),-5*x(4)];
gh3 = [3*x(1)^2,3*x(2)^2,0,0,0];
Dh = [gh1;gh2;gh3];
gh = Dh';