clear;clc
Path = pwd;
addpath(Path);
cd Ex3.10
%% Three stationary points
x0 = [0; 0; 0; 0];
[x_star,Fx] = FMSQP(@funf,@gradf,@func,@gradc,x0)

x0 = [0; -1; 1; -1];
[x_a,Fx] = FMSQP(@funf,@gradf,@func,@gradc,x0)

x0 = [ 0; -1; 1; 1];
[x_b,Fx] = FMSQP(@funf,@gradf,@func,@gradc,x0)


%% Output  for Table 2
opts.varbose = 1;
opts.nmax = 9;
opts.epsilon = 1e-10;

fprintf('Output for x_0^1\n\n')
x01 = [5; 5; 5; 5];
[x,Fx] = FMSQP(@funf,@gradf,@func,@gradc,x01,opts)

fprintf('\n')
fprintf('Output for x_0^2\n')
x02 = [-5; -5; -5; -5];
[x,Fx] = FMSQP(@funf,@gradf,@func,@gradc,x01,opts)

%% Print for Fig.1
run('PrintFig.m')
cd(Path)
rmpath(Path);
