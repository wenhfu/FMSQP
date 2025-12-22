clear;clc
Path = pwd;
addpath(Path);
cd Ex5.1
fid = fopen("Output_Loop.txt",'a');
opts.epsilon = 1e-4;
for t = [1,5] % [1,5,10,50,100]
    Fx_multi = [];
    Fx_single = [];
    n = 2;
    Nmax = 10; % 1000
    toc_multi = 0;
    nit_multi = zeros(Nmax,1);
    nf_multi = zeros(Nmax,1);
    ng_multi = zeros(Nmax,1);
    toc_single = 0;
    nit_single = zeros(Nmax,1);
    nf_single = zeros(Nmax,1);
    ng_single = zeros(Nmax,1);
    for Nit = 1:Nmax
        L = 2;
        rng('shuffle');
        X0 = t*rand(n,L)-2*t;
        % [~,Fx,output] = FMSQP(X0,opts);
        [~,Fx,output] = FMSQP(@funf,@gradf,@func,@gradc,X0,opts);
        toc_multi = toc_multi + output.time;
        nit_multi(Nit) = output.nit;
        nf_multi(Nit) = output.nf;
        ng_multi(Nit) = output.ng;
        if ~isempty(Fx)
            Fx_multi = [Fx_multi;Fx];
        end
        I_multi = find(abs(Fx_multi-1)<1e-3);
        Pre_multi = length(I_multi)/length(Fx_multi)*100;
        time_multi = toc_multi/Nit;
        L = 1;
        rng('shuffle');
        X0 = t*rand(n,L)-2*t;
        % [~,Fx,output] = FMSQP(X0,opts);
        [~,Fx,output] = FMSQP(@funf,@gradf,@func,@gradc,X0,opts);
        toc_single = toc_single + output.time;
        nit_single(Nit) = output.nit;
        nf_single(Nit) = output.nf;
        ng_single(Nit) = output.ng;
        if ~isempty(Fx)
            Fx_single = [Fx_single;Fx];
        end
        I_single = find(abs(Fx_single-1)<1e-3);
        Pre_single = length(I_single)/length(Fx_single)*100;
        time_single = toc_single/Nit;
        fprintf('\nt = %4d, Nit =%4d: \n',t, Nit);
        fprintf('Pre_multi=%.2f%%, time_multi=%.2f \n',Pre_multi,time_multi);
        fprintf('Pre_single=%.2f%% time_single=%.2f \n',Pre_single,time_single);
    end
    % Print for Table 3
    fprintf(fid,'\\multirow[c]{2}{*}{$x_0\\in U(-%d,%d)$}\n',t,t);
    fprintf(fid,'& 1 & %4.2f & %4.2f & %4.2f & %4.2f & %4.2f\\%% \\\\\n',mean(nit_single),mean(nf_single),mean(ng_single),toc_single/Nmax,Pre_single);
    fprintf(fid,'& 2 & %4.2f & %4.2f & %4.2f & %4.2f & %4.2f\\%% \\\\\n',mean(nit_multi),mean(nf_multi),mean(ng_multi),toc_multi/Nmax,Pre_multi);
end
cd(Path)








