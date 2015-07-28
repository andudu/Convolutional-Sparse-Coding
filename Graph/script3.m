
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Load  Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load a saved noise
%load('noise_data.mat');
%snoise = s;
% snoise = randn(256,256)*.1;
s_ref = single(stdimage('lena.grey')) / 255;
% snoise = randn(size(s_ref))*.1;
% s = s_ref+snoise;
wsz = [60,60];
psz = [12,12];
% neig = 30;
% [sl,sh] = lowpass(s,5,15);

load('lena.mat'); %because I have no license for some reason.
[L,sh] = graphgen(sh,wsz,psz,neig);
disp('graph generated');
sl = sl(1:size(sh,1),1:size(sh,2));
s_ref = s_ref(1:size(sh,1),1:size(sh,2));


% Load dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');



%%%%%%%%%%%%%%%%%%%%%%% set up parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

 mu_all = 0.1:0.1:1;
%mu_all = 0.1:0.1:0.2;
opt = {};
opt.Verbose = 0;
opt.MaxMainIter = 300;
opt.rho =10;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;

%optimize
lambda_all = 0.1:0.02:0.3;
%lambda_all = 0.2:0.01:0.21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% optimizing over all param %%%%%%%%%%%%%%%%%%%%
% reconstructing conv
scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

pcn_rec = zeros(size(lambda_all));
pnl_rec = zeros(length(lambda_all), length(mu_all));
pcn_max = 0;
pnl_max = 0;
lambdacn_max = 0;
lambdanl_max = 0;
munl_max = 0;



for i = 1:length(lambda_all)
    disp([num2str(i),' out of ',num2str(length(lambda_all))]);
    lambda = lambda_all(i);
    [Xcn,~] = cbpdn(D,sh,lambda,opt);
    shreccn = scnv(D,Xcn);
    sreccn = shreccn+sl;
    p_cn = psnr(sreccn,s_ref);
    pcn_rec(i) = p_cn;
    if p_cn > pcn_max,
        pcn_max = p_cn;
        lambdacn_max = lambda;
    end
    
    for j = 1:length(mu_all)
        mu = mu_all(j);
        [Xnl,~]= cbpdn_L(D,sh,L,lambda,mu,opt);
        shrecnl = scnv(D,Xnl);
        srecnl = shrecnl+sl;
        p_nl = psnr(srecnl,s_ref);
        pnl_rec(i,j) = p_nl;
        if p_nl > pnl_max
            pnl_max = p_nl;
            lambdanl_max = lambda;
            munl_max = mu;
        end
    end
end

tag = 'parmsearchmu_12_512.mat';
save(['ParamNL/',tag],'pcn_rec','pnl_rec','pcn_max','pnl_max','lambdacn_max','lambdanl_max','munl_max','mu_all','lambda_all');