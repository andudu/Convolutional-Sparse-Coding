%this script compares sparse coding between applying the filter 
%and ways to avoid it. 

% Load a standard dictionary
load([sporco_path '/Data/ConvDict.mat']);
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
num_dict = size(D,3);


% Load test image
s = single(rgbtogrey(stdimage('barbara')))/255; %use barbara here
if isempty(s),
  error('Data required for demo scripts has not been installed.');
end

scnv = @(d,x) ifft2(sum(bsxfun(@times, fft2(d, size(x,1), size(x,2)), ...
                               fft2(x)),3), 'symmetric');

% Highpass filter test image
npd = 16;
fltlmbd = 5;
[sl, sh] = lowpass(s, fltlmbd, npd);


%setting up parameters
lambda = 0.05;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 500;
opt.rho = 100*lambda + 1;
opt.RelStopTol = 1e-3;
opt.AuxVarObj = 0;
opt.HighMemSolve = 1;



% Compute representation using cbpdngr with no pass but high regularizer
mu = 1;
opt.L1Weight = reshape([ones(1,num_dict-1),.8],1,1,num_dict);
opt.GrdWeight = reshape([zeros(1,num_dict-1),10],1,1,num_dict);
[X2, optinf2] = cbpdngr(D, s, lambda, mu, opt);


% Compute representation using cbpdn
opt = rmfield(opt,'GrdWeight');
opt.L1Weight = 1;
[X1, optinf1] = cbpdn(D, sh, lambda, opt);


%save in sparse format
X1 = compress_nd(X1);
X2 = compress_nd(X2);
save(strcat(sporco_path,'/Results/benchmark.mat'),'X1','optinf1');
save(strcat(sporco_path,'/Results/grdx.mat'),'X2','optinf2');






