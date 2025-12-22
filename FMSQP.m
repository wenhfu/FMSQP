function [x,fx,output] = FMSQP(X0,opts)
% [x,fval,exitflag,output,lambda] = NPFopt(funf,func,x0,MF,opts)
% function [x,fval,exitflag,output,lambda] = FMSQP(fx,cx,X0,opts)
gamma = 0.6; eps_u = 10; eta = 1e-3; epsilon = 1e-5; 
nmax = 500; sigma0 = 1; sigma = sigma0; 
nit = 0; nf = 0; ng = 0; varbose = 0;
if nargin == 2
% if nargin == 4
    if isfield(opts, 'varbose'); varbose = opts.varbose; end
    if isfield(opts, 'epsilon'); epsilon = opts.epsilon; end
    if isfield(opts, 'nmax'); nmax = opts.nmax; end
end

tic

% [funf,gradf] = @(x)fx(x)
% [func,gradc] = @(x)cx(x)

[n,L] = size(X0);
for j=1:L
    xold = X0(:,j);
    Xold{j} = xold;
    Fxold(j) = funf(xold); Gfxold{j} = gradf(xold); nf = nf + 1; ng = ng + 1;
    Cxold{j} = func(xold); Gcxold{j} = gradc(xold);
    Vxold(j) = norm(max(Cxold{j},0),'inf');
    Bk{j} = eye(n); Hk{j} = eye(n); Sigma(j) = sigma;
end
epsM = 1e-4; beta = 20;
Jmax = floor(log(max(Vxold+0.1)/epsM)/log(beta));
M = beta.^(1:Jmax)*epsM;
M = [0,M];

J_succeed = [];
J_fail = [];
J_old = 1:L;
Vcut = [];
Data_table = [];
while nit < nmax
    Jk_opt = [];
    Vcut = Vxold;
    for i=1:length(M)-1
        I1 = find(Vcut>=M(i));
        I2 = find(Vcut<=M(i+1));
        Iv = setdiff(intersect(I1,I2),J_fail);
        if ~isempty(Iv) && length(Iv) >= 2 && nit > 0
            F_ki = min(Fxold(Iv));
            If = find(Fxold>F_ki);
            I = intersect(Iv,If);
            Jk_opt = union(Jk_opt,I);
            J_fail = union(J_fail,I);
        end
    end
    for j = setdiff(J_old,Jk_opt)
        fxold = Fxold(j); gfxold = Gfxold{j};
        cxold = Cxold{j}; gcxold = Gcxold{j};
        vxold = Vxold(j);
        shift = sub1(eye(n)/max(100,vxold^2),cxold,gcxold);
        [dk,lambda,shift] = sub2(gfxold,Bk{j},cxold,gcxold,shift);
        Norm_dk(j) = norm(dk);
        Lkv{j,nit+1} = shift.lkv_dkopt;
        Dk{j,nit+1} = dk;
        Vk(j,nit+1) = vxold;
        Lambda{j} = lambda;
        if norm(dk,'inf') <= epsilon
            J_succeed = union(J_succeed,j);
        elseif isempty(dk)
            J_fail = union(J_fail,j);
        else
            uk = shift.lkv_dkopt/norm(dk);
            Uk(j,nit+1) = uk;
        end
        if varbose == 1
            fprintf('%2d & %.4f & %.4f & %.1e & %.1e \n',nit,fxold,vxold,norm(dk,'inf'),uk);
        end
    end
    if isempty(setdiff(1:L,union(J_succeed,J_fail)))
        break
    end
    hat_u = mean(Uk(:,nit+1));
    hat_v = mean(Vk(:,nit+1));
    Jk_fea = [];
    for j = setdiff(J_old,union(Jk_opt,J_succeed))
        if (Uk(j,nit+1) > hat_u + eps_u) && (Vk(j,nit+1) > hat_v) && nit > 0% && length(setdiff(J_old,union(Jk_opt,J_succeed)))==1
            if length(setdiff(J_old,union(J_fail,J_succeed))) > 1
                Jk_fea = union(Jk_fea,j);
                J_fail = union(J_fail,j);
                continue
            end
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
        if varbose == 2
            fprintf('%2d & %2d & %.2e & %.2f & %.4e & %.4e & %.4f\n',nit,j,norm(dk,'inf'),alpha,shift.lkv_dkopt,vxnew,fxnew);
        end

        % lkv = shift.lkv_dkopt;
        % uxold = lkv/norm(dk);
        % norm(uxold - Uk(j,nit+1)) % 不相等
        % Data_table(nit+1,:) = [nit,fxold,vxold,norm(dk),lkv,uxold]
        % Data_table(nit+1,:) = [nit,fxold,vxold,norm(dk),lkv,Uk(j,nit+1)]
        % fprintf('%d & %.4f & %.e & %.0e & %.e & %.0e \\\\\n',nit,fxold,vxold,norm(dk),lkv,uxold);
        if alpha <= 1e-3
            if abs(fxold-fxnew)/abs(fxnew)<1e-5 && vxnew < 1e-5
                J_succeed = union(J_succeed,j);
            else
                if length(setdiff(J_old,union(J_fail,J_succeed))) > 1
                    J_fail = union(J_fail,j);
                    continue
                end
            end
        end
        gfxnew = gradf(xnew);
        gcxnew = gradc(xnew); ng = ng + 1;
        s = xnew - xold;
        lambda = Lambda{j};
        q = gfxnew + gcxnew*lambda - ( gfxold + gcxold*lambda);
        Bk{j} = bfgs(Bk{j},s,q);
        mm = eig(Bk{j});
        if cond(Bk{j})>1e+6; Bk{j}=eye(n); end
        if min(mm)<1e-6; Bk{j} = Bk{j}+eye(n); end
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
m = length(Cxold{1});
J_feasible = find(Vxold(J_succeed)<1e-3);
[fx,J] = min(Fxold(J_succeed(J_feasible)));
if isempty(J)
    vx = min(Vxold);
    J_feasible = find(Vxold == vx);
    J = J_feasible(1);
    fx = Fxold(J);
    x = Xold{J};
else
    vx = Vxold(J_succeed(J_feasible(J)));
    x = Xold{J};
end
output.m = m;
output.vx = vx;
output.nit = nit;
output.nf = nf;
output.ng = ng;
output.time = toc;
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
[Sol,~,exitflag] = solve(prob,Sol0,options=Options);
shift.dkfea = Sol.d;
shift.cons_QP12 = cxold+gcxold'*shift.dkfea;
shift.lkv_dkfea = norm(max(shift.cons_QP12,0),'inf');
shift.rkfea = Sol.r;
shift.flag = 0;
if exitflag<=0
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
