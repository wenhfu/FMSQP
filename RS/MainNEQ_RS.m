clear;  clc
t = 5;
X0 = [t*ones(1,4);-t*ones(1,4)]';
n = 4;
gamma = 0.6; eta = 1e-3; eps_u = 10;
nmax = 500; epsilon = 1e-5; sigma0 = 1; sigma = sigma0;
nit = 0; flag = 0; nf = 0; ng = 0;

L = size(X0,2);
% L = 10;
% rng(1);
% rng('shuffle');
% X0 = rand(n,L);
% X0 = 2*randn(n,L);
% X0 = ones(n) - 2*eye(n); L = n;
% X0 = ones(n) - eye(n); L = n;
% X0 = ones(n) + eye(n); L = n;
% X0 = ones(n) - 2*eye(n).*diag(randn(n,1)); L = n;
Hess_c1 = 2*eye(4);
Hess_c2 = -2*eye(4);
Hess_c3 = diag([2,4,2,4]);
Hess_c4 = -diag([2,4,2,4]);
Hess_c5 = diag([4,2,2,0]);
Hess_c6 = -diag([4,2,2,0]);
% Hess_c = diag([2,2,2,2,-2,-2,-2,-2,2,4,2,4,-2,-4,-2,-4,4,2,2,0,-4,2,2,0]);
for j=1:L
    xold = X0(:,j);
    Xold{j} = xold;
    Fxold(j) = funf(xold); Gfxold{j} = gradf(xold); nf = nf + 1; ng = ng + 1;
    Cxold{j} = func(xold); Gcxold{j} = gradc(xold);
    Vxold(j) = norm(max(Cxold{j},0),'inf');
    Bk{j} = eye(n); Sigma(j) = sigma;
    Hk{j} = eye(n);
end
eps_Uk = 0.01; epsM = 1e-4; beta = 20;
Jmax = floor(log(max(Vxold+0.1)/epsM)/log(beta));
M = beta.^(1:Jmax)*epsM;
M = [0,M];

J_succeed = [];
J_fail = [];
J_old = 1:L;
Vcut = [];
idX = [];
while 1
    Jk_opt = [];
    for j = setdiff(J_old,Jk_opt)
        fxold = Fxold(j); gfxold = Gfxold{j};
        cxold = Cxold{j}; gcxold = Gcxold{j};
        vxold = Vxold(j);
        shift = sub1(Hk{j},cxold,gcxold);
        [dk,lambda,shift] = sub2(gfxold,Bk{j},cxold,gcxold,shift);
        Hk{j} = Hess_c1*shift.lambda(1)+Hess_c2*shift.lambda(2)+Hess_c3*shift.lambda(3)+Hess_c4*shift.lambda(4)+Hess_c5*shift.lambda(5)+Hess_c6*shift.lambda(6);
        mm = min(eig(Hk{j}));
        while mm < 1e-4
            Hk{j} = Hk{j}+eye(n);
            mm = min(eig(Hk{j}));
        end
        Norm_dk(j) = norm(dk);
        Lkv{j,nit+1} = shift.lkv_dkopt;
        Dk{j,nit+1} = dk;
        Vk(j,nit+1) = vxold;
        Fk(j,nit+1) = fxold;
        Lambda{j} = lambda;
        Uk(j,nit+1) = shift.lkv_dkopt/norm(dk);
        if norm(dk,'inf') <= epsilon
            J_succeed = union(J_succeed,j);
            idX(j) = nit+1;
        % elseif isempty(dk)
        %     J_fail = union(J_fail,j);
        % else
            % Uk(j,nit+1) = shift.lkv_dkopt/norm(dk);
        end
    end
    if isempty(setdiff(1:L,union(J_succeed,J_fail)))
        break
    end
    hat_u = mean(Uk(:,nit+1));
    hat_v = mean(Vk(:,nit+1));
    Jk_fea = [];
    for j = J_old%setdiff(J_old,union(Jk_opt,J_succeed))
        if (Uk(j,nit+1) > hat_u + eps_u) && (Vk(j,nit+1) > hat_v) && nit > 0 % && length(setdiff(J_old,union(Jk_opt,J_succeed)))==1
            % if length(setdiff(J_old,union(J_fail,J_succeed))) > 1
                % Jk_fea = union(Jk_fea,j);
                % J_fail = union(J_fail,j);
                fprintf('CUT:%d \n',j);
                % idX(j) = nit;
                % continue
            % end
        end
        xold = Xold{j}; dk = Dk{j,nit+1};
        fxold = Fxold(j); gfxold = Gfxold{j};
        cxold = Cxold{j}; gcxold = Gcxold{j};
        vxold = Vxold(j); sigma = Sigma(j);
        alpha = 1;
        sigma = max(sigma,norm(Lambda{j},1)+1);
        Sigma(j) = sigma;
        Pxold = fxold + sigma*vxold;
        while alpha > 1e-6
            xnew = xold + alpha*dk;
            fxnew = funf(xnew); nf = nf + 1;
            cxnew = func(xnew);
            vxnew = norm(max(cxnew,0),'inf');
            Pxnew = fxnew + sigma*vxnew;
            eqn = Pxnew - Pxold + eta*alpha*(dk'*Bk{j}*dk+(Lkv{j}-Vk(j,nit+1))); 
            if eqn <= 0
                break
            end
            alpha = alpha * gamma;
        end
        fprintf('%2d & %2d & %.2e & %.2f & %.4e & %.4e & %.4f\n',nit,j,norm(dk,'inf'),alpha,shift.lkv_dkopt,vxnew,fxnew);
        gfxnew = gradf(xnew);
        gcxnew = gradc(xnew);ng = ng + 1;
        s = xnew - xold;
        lambda = Lambda{j};
        q = gfxnew + gcxnew*lambda - ( gfxold + gcxold*lambda);
        Bk{j} = bfgs(Bk{j},s,q);
        mm = eig(Bk{j});
        if cond(Bk{j})>1e+6
            Bk{j}=eye(n);
        end
        if min(mm)<1e-6
            Bk{j} = Bk{j}+eye(n);
        end
        % q = gcxnew*lambda - ( gcxold*lambda);
        % Hk{j} = bfgs(Hk{j},s,q);
        % mm = eig(Hk{j});
        % if cond(Hk{j})>1e+6
        %     Hk{j}=eye(n);
        % end
        % if min(mm)<1e-6
        %     Hk{j} = Hk{j}+eye(n);
        % end
        xold = xnew;
        fxold = fxnew; gfxold = gfxnew;
        cxold = cxnew; gcxold = gcxnew;
        vxold = vxnew;
        Xold{j} = xold;
        Fxold(j) = fxold; Gfxold{j} = gfxold;
        Cxold{j} = cxold; Gcxold{j} = gcxold;
        Vxold(j) = vxold;
    end
    J_old = setdiff(J_old,union(J_fail,J_succeed));
    nit = nit + 1;
end
j0 = find(Vxold(J_succeed)==min(Vxold(J_succeed)));
m = length(Cxold{1});
J_feasible = find(Vxold(J_succeed)<1e-3);
[f_num,J] = min(Fxold(J_succeed(J_feasible)));
v_num = Vxold(J_succeed(J_feasible(J)));
Fx = fxold;
x = xold;

% figure
% hold on
% Len = 7;
% U{1} = Uk(1,max(1,idX(1)-Len+1):idX(1));
% U{2} = Uk(2,max(1,idX(2)-Len+1):idX(2));
% U{i} = Uk(i,1:Len);
% for i = 1:L
%     plot((U{i}));
% end
% axis([1,Len,-1,1000])
fprintf('\n');
for t = 1:max(idX)
    % if Vk(1,t)>1e-5
        fprintf('%1d & %.4f & %.4f & %.1e & %.1e && %.4f & %.4f & %.1e & %.1e \\\\\n',t-1,Fk(1,t),Vk(1,t),norm(Dk{1,t},'inf'),Uk(1,t),Fk(2,t),Vk(2,t),norm(Dk{2,t},'inf'),Uk(2,t));
    % else
        % fprintf('%1d & - & - & - & - & %.2f & %.2f & %.2e & %.2f \\\\\n',t-1,Fk(2,t),Vk(2,t),norm(Dk{2,t},'inf'),Uk(2,t));
    % end
end














%% Subproblems

function B=bfgs(B,d,y)
bd=B*d;dbd=d'*bd;dy=d'*y;
if (dy<0.2*dbd)
    theta=0.8*dbd/(dbd-dy);
    y=theta*y+(1.0-theta)*bd;
end
dybar=d'*y;
B=B+y*y'/dybar-bd*bd'/dbd;
end

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
[Sol,~,exitflag,~,Lambda] = solve(prob,Sol0,options=Options);
shift.dkfea = Sol.d;
shift.cons_QP12 = cxold+gcxold'*shift.dkfea;
shift.lkv_dkfea = norm(max(shift.cons_QP12,0),'inf');
shift.rkfea = Sol.r;
shift.flag = 0;
shift.lambda = Lambda.Constraints.c3;
if exitflag<=0
    % disp('sub1 fails')
    % disp(exitflag)
    shift.flag = -1;
end
end

function [dk,lambda,shift] = sub2(gfxold,Bk,cxold,gcxold,shift)
[n,mi] = size(gcxold);
rkfea = shift.rkfea;
prob = optimproblem;
d = optimvar('d',n,1);
Sol0.d = zeros(n,1);
prob.Objective = gfxold'*d+0.5*d'*Bk*d;
prob.Constraints.c3 = cxold+gcxold'*d<=rkfea;
Options = optimset('Display','none');
[Sol,~,exitflag,~,Lambda] = solve(prob,Sol0,options=Options,solver='quadprog');
shift.dkopt = Sol.d;
if exitflag<=0 || isempty(Lambda.Constraints.c3)
    % disp('sub2 fails')
    % disp(exitflag)
    shift.flag = -2;
    dk = [];
    lambda = [];
    shift.lkv_dkopt = [];
else
    dk = Sol.d;
    y = Lambda.Constraints.c3;
    lambda = y;
    shift.cons_QP22 = cxold+gcxold'*dk;
    shift.lkv_dkopt = norm(max(shift.cons_QP22,0),'inf');
    shift.dkopt = Sol.d;
end
end
