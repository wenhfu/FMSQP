clear; warning('off');  clc
Path = pwd;
addpath(Path);
Example{1}.Name = "G1"; Example{1}.n = 13; Example{1}.fstar = -15.00;
Example{2}.Name = "G2"; Example{2}.n = 20; Example{2}.fstar = -0.8036; %***** not an NLP
Example{3}.Name = "G3"; Example{3}.n = 10; Example{3}.fstar = -1.00;
Example{4}.Name = "G4"; Example{4}.n = 5; Example{4}.fstar = -30665.539;
Example{5}.Name = "G5"; Example{5}.n = 4; Example{5}.fstar = 5126.498;
Example{6}.Name = "G6"; Example{6}.n = 2; Example{6}.fstar = -6961.814;
Example{7}.Name = "G7"; Example{7}.n = 10; Example{7}.fstar = 24.306;
Example{8}.Name = "G8"; Example{8}.n = 2; Example{8}.fstar = -0.095825;
Example{9}.Name = "G9"; Example{9}.n = 7; Example{9}.fstar = 680.63;
Example{10}.Name = "G10"; Example{10}.n = 8; Example{10}.fstar = 7049.331; % ***** failed 
Example{11}.Name = "G11"; Example{11}.n = 2; Example{11}.fstar = 0.75;
Example{12}.Name = "G12"; Example{12}.n = 3; Example{12}.fstar = -1.00; % not an NLP
Example{13}.Name = "G13"; Example{13}.n = 5; Example{13}.fstar = 0.05395;

opts.varbose = 0;
Fail_idx = [2,12];
for idx = 1:13
    fprintf("Name: %s                  Best value: %12.4f\n",Example{idx}.Name,Example{idx}.fstar);
    n = Example{idx}.n;
    rng(1); L = 10;
    X0 = rand(n,L);
    Path = fullfile(pwd,Example{idx}.Name);
    addpath(Path);
    if ismember(idx,Fail_idx)
        fprintf('  not an NLP\n')
        Example{idx}.m = [];
        Example{idx}.f_num = [];
        Example{idx}.v_num = [];
        Example{idx}.nit = [];
        Example{idx}.nf = [];
        Example{idx}.ng= [];
    else
        [x,fx,output] = FMSQP(@funf,@gradf,@func,@gradc,X0,opts);
        Example{idx}.f_num = fx;
        Example{idx}.x = x;
        Example{idx}.m = output.m;
        Example{idx}.v_num = output.vx;
        Example{idx}.nit = output.nit;
        Example{idx}.nf = output.nf;
        Example{idx}.ng= output.ng;
    end
    rmpath(Path);
    fprintf('\n')
end

save Data_G Example Fail_idx

f_OSCARS = {'-14.99946';'-0.628909';'-0.999611';'-30655.66';'5043';'-6962.046';'24.31826';'-0.095825';'680.6311';'7049.638';'0.7499';'-1';'0.0539428';'-47.75237';'961.71502';'-1.9051686';'-0.865944'};
% for idx = 1:13
%     if ~ismember(idx,Fail_idx)
%         fprintf('%9s &%4d &%4d &%4d &%4d &%4d &%12.4e &%12.4f &%12.4f &%12s \\\\ \n',Example{idx}.Name,Example{idx}.n,Example{idx}.m,Example{idx}.nit,Example{idx}.nf,Example{idx}.ng,Example{idx}.v_num,Example{idx}.f_num,Example{idx}.fstar,f_OSCARS{idx});
%     end
% end

for idx = 1:13
    if ~ismember(idx,Fail_idx)
        fprintf('%9s &%4d &%4d &%4d &%12.4e &%12.4f &%12.4f &%12s \\\\ \n',Example{idx}.Name,Example{idx}.nit,Example{idx}.nf,Example{idx}.ng,Example{idx}.v_num,Example{idx}.f_num,Example{idx}.fstar,f_OSCARS{idx});
    end
end

