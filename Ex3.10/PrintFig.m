% clear;clc
n = 4;
EPS0 = [1e-0,1e-1,1e-2,1e-3,1e-4,1e-5];
Data_u = [];

fprintf('x_star:\n')
x_star = [0,1,2,-1]'; % x_star
k=1;
for eps0 = EPS0
xold = x_star+eps0;
fxold = funf(xold); gfxold = gradf(xold);
cxold = func(xold); gcxold = gradc(xold);
vxold = norm(max(cxold,0),'inf');
shift = sub1(eye(n),cxold,gcxold);
[dk,lkv] = subQP(gfxold,eye(n),cxold,gcxold,shift);
uxold = lkv/norm(dk);
Delta_u(k,1) = uxold; k=k+1;
% fprintf('eps0=%.2e, fk=%.2e, vk=%.2e, dk=%.2e, lkv=%.2e, lkv/dk=%.2e\n',eps0,fxold,vxold,norm(dk),lkv,uxold);
fprintf('%.0e & %.4f & %.e & %.0e & %.e & %.0e \\\\\n',eps0,fxold,vxold,norm(dk),lkv,uxold);
end

fprintf('\nx_a:\n')
x_a = [-0.0159,-1.1277,1.1118,-1.3800]'; % x_a
k=1;
for eps0 = EPS0
xold = x_a+eps0;
fxold = funf(xold); gfxold = gradf(xold);
cxold = func(xold); gcxold = gradc(xold);
vxold = norm(max(cxold,0),'inf');
shift = sub1(eye(n),cxold,gcxold);
[dk,lkv] = subQP(gfxold,eye(n),cxold,gcxold,shift);
uxold = lkv/norm(dk);
Delta_u(k,2) = uxold; k=k+1;
% fprintf('eps0=%.2e, fk=%.2e, vk=%.2e, dk=%.2e, lkv=%.2e, lkv/dk=%.2e\n',eps0,fxold,vxold,norm(dk),lkv,uxold);
fprintf('%.0e & %.4f & %.4f & %.0e & %.4f & %.0e \\\\\n',eps0,fxold,vxold,norm(dk),lkv,uxold);
end

fprintf('\nx_b:\n')
x_b = [-0.0351,-1.5525,1.5174,1.2325]'; % x_b
k=1;
for eps0 = EPS0
xold = x_b+eps0;
fxold = funf(xold); gfxold = gradf(xold);
cxold = func(xold); gcxold = gradc(xold);
vxold = norm(max(cxold,0),'inf');
shift = sub1(eye(n),cxold,gcxold);
[dk,lkv] = subQP(gfxold,eye(n),cxold,gcxold,shift);
uxold = lkv/norm(dk);
Delta_u(k,3) = uxold; k=k+1;
% fprintf('eps0=%.2e, fk=%.2e, vk=%.2e, dk=%.2e, lkv=%.2e, lkv/dk=%.2e\n',eps0,fxold,vxold,norm(dk),lkv,uxold);
fprintf('%.0e & %.4f & %.4f & %.0e & %.4f & %.0e \\\\\n',eps0,fxold,vxold,norm(dk),lkv,uxold);
end

fig = figure;
fig.Position = [600,300,500,300];
grid on
hold on
plot((Delta_u(:,2)),'ro-.','LineWidth',1.5)
plot((Delta_u(:,3)),'bsquare--','LineWidth',1.5)
plot((Delta_u(:,1)),'k*:','LineWidth',1.5)
axis([1,6,-50,500])
xlabel('$\varepsilon$','Interpreter','latex');
xticklabels({'10^{0}','10^{-1}','10^{-2}','10^{-3}','10^{-4}','10^{-5}'})
ylabel('$u(x,d)$','Interpreter','latex');
legend({'$x_a+\varepsilon$','$x_b+\varepsilon$','$x_*+\varepsilon$'},'Interpreter','latex','Location','northwest')
print('fig3_1', '-dpng', '-r600')


%% Subproblems

function shift = sub1(Bk,cxold,gcxold)
[n,mi] = size(gcxold);
prob = optimproblem;
d = optimvar('d',n,1);
r = optimvar('r',1,1);
Sol0.d = zeros(n,1);
Sol0.r = zeros(1,1);
prob.Objective = r+0.5*d'*Bk*d;
prob.Constraints.c3 = cxold+gcxold'*d<=r;
prob.Constraints.c4 = r>=0;
Options = optimset('Display','none');
Sol = solve(prob,Sol0,options=Options);
shift.rkfea = Sol.r;
end

function [dk,lkv] = subQP(gfxold,Bk,cxold,gcxold,shift)
[n,~] = size(gcxold);
prob = optimproblem;
d = optimvar('d',n,1);
Sol0.d = zeros(n,1);
prob.Objective = gfxold'*d+0.5*d'*Bk*d;
prob.Constraints.c1 = cxold+gcxold'*d<=shift.rkfea;
Options = optimset('Display','none');
Sol = solve(prob,Sol0,options=Options,solver='quadprog');
dk = Sol.d;
lkv = norm(max(cxold+gcxold'*dk,0),'inf');
end
