clear;clc; warning('off')
Path = pwd;
addpath(Path);
addpath(fullfile(Path,'/ExCUTE'));
addpath('~/.local/matcutest/mtools/src');


pname_hs = table2cell(readtable('PNAME_HS_ALL')); N = size(pname_hs,1);
Fail_idx = [];
Data_FMSQP_1 = cell(length(pname_hs),9);
Data_FMSQP_2 = cell(length(pname_hs),9);
Data_NaturalSQP = cell(length(pname_hs),9);
Data_LSBO = cell(length(pname_hs),9);
Data_SQuID = cell(length(pname_hs),9);
Data_SteeringSQP = cell(length(pname_hs),9);
Data_ALMSQP = cell(length(pname_hs),9);
opts.varbose = 0;
% N_time = 1;
% for name_idx = 17
for name_idx = 1 : length(pname_hs)
% for name_idx = [68,75,77,100,111,112,118,119]
% for name_idx = [77,78]
    pname = pname_hs{name_idx};
    prob = macup(pname);
    x0 = prob.x0;
    fx = @(x) funf(x,prob);
    gfx = @(x) gradf(x,prob);
    cx = @(x) func(x,prob);
    gcx = @(x) gradc(x,prob);
    fprintf('\nname_idx = %3d:  Name = %9s,   \n',name_idx,pname);
    %% FMSQP_SigleStart
    [~,~,output] = FMSQP(fx,gfx,cx,gcx,x0,opts); Data_FMSQP_1(name_idx,:) = [name_idx,pname,{output.exitflag},{output.nit},{output.nf},{output.ng},{output.fx},{output.vx},{output.time}];
    %% FMSQP_MultiStart
    [~,~,output] = FMSQP(fx,gfx,cx,gcx,[x0,(1e-8)*randn(size(x0)),ones(size(x0)),-ones(size(x0))],opts); Data_FMSQP_2(name_idx,:) = [name_idx,pname,{output.exitflag},{output.nit},{output.nf},{output.ng},{output.fx},{output.vx},{output.time}];
end
% save Data_FMSQP_1 Data_FMSQP_1
% save Data_FMSQP_2 Data_FMSQP_2
I_single = find(cell2mat(Data_FMSQP_1(:,3))<2); Data_FMSQP_1(I_single,:) % 8
I_multi = find(cell2mat(Data_FMSQP_2(:,3))<2); Data_FMSQP_2(I_multi,:) % 7
I_success = [];
for i=1:size(Data_FMSQP_2,1)
    if Data_FMSQP_2{i,3} == 2 && Data_FMSQP_1{i,8} < 1e-3 && Data_FMSQP_2{i,8} < 1e-3
        if abs(Data_FMSQP_1{i,7}-Data_FMSQP_2{i,7}) / max(1,abs(Data_FMSQP_1{i,7})) > 0.01 && Data_FMSQP_1{i,7} >= 0.99*Data_FMSQP_2{i,7}
            I_success = union(I_success,i);
        end
    end
end
I_success

I = [I_success,setdiff(I_single,I_multi)];
for i = I
    fprintf('\\midrule \\multirow{2}{*}{%9s} \n',Data_FMSQP_1{i,2})
    fprintf('& 1  & %4d & %4d & %4d & %9.4f & %9.4e \\\\ \n',Data_FMSQP_1{i,4:8})
    fprintf('& 10 & %4d & %4d & %4d & %9.4f & %9.4e \\\\ \n',Data_FMSQP_2{i,4:8})
end



%% sub_functions
function f = funf(x,prob)
f = prob.objective(x);
end

function gf = gradf(x,prob)
[~,gf] = prob.objective(x);
end

function c = func(x,prob)
n = length(x);
if prob.numnlcon > 0
    [g,h] = prob.nonlcon(x);
else
    g = []; h = [];
end
if prob.numlineq > 0
    g = [g;prob.Aineq*x-prob.bineq];
end
if prob.numleq > 0
    h = [h;prob.Aeq*x-prob.beq];
end
for i=1:n
    if ~isempty(prob.lb)
        if ~isinf(prob.lb(i))
            g = [g;-x(i)+prob.lb(i)];
        end
        if ~isinf(prob.ub(i))
            g = [g;x(i)-prob.ub(i)];
        end
    end
end
c = [h;-h;g];
end

function gc = gradc(x,prob)
n = length(x);
if prob.numnlcon > 0
    [~,~,Jg,Jh] = prob.nonlcon(x);
else
    Jg = []; Jh = [];
end
if prob.numlineq > 0
    Jg = [Jg,prob.Aineq'];
end
if prob.numleq > 0
    Jh = [Jh,prob.Aeq'];
end
for i=1:n
    if ~isempty(prob.lb)
        if ~isinf(prob.lb(i))
            Jgi = zeros(n,1); Jgi(i) = -1;
            Jg = [Jg,Jgi];
        end
        if ~isinf(prob.ub(i))
            Jgi = zeros(n,1); Jgi(i) = 1;
            Jg = [Jg,Jgi];
        end
    end
end
gc = [Jh,-Jh,Jg];
end