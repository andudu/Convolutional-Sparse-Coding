% script for testing the fake image inpainting problem
% lambda = .01 for the reconstruction. 
% cbpdn is used for both the Laplacian and Convolutional Dictionaries


perc_noise = .5;
lambda_inp = .02;

% load the image peppers
S = double(rgb2gray(stdimage('peppers')))/255;
maxiter = 300;

num_i = 3;
num_j = 3;


mu_all = [];
lambda_all = [];
num1 = [];
num2 = [];
psnr_rec = zeros(num_i,num_j,2);
snr_rec = zeros(num_i,num_j,2);

for i = 1:num_i %this needs some adjusting
    for j = 1:num_j      
        disp([num2str(i),',',num2str(j)]);
        %load the dictionary
        load(['DictComp',num2str(i),num2str(j)]);
        mu_all(i,j) = mu;
        lambda_all(i,j) = lambda;
        
        %prune the dictionaries        
        Aind1 = [];
        for k = 1: size(D1,3)
            if(nnz(X1(:,:,i,:))>=20)
                Aind1 = [Aind1,i];
            end
        end        
        Aind2 = [];
        for k = 1: size(D2,3)
            if(nnz(X2(:,:,i,:))>=20)
                Aind2 = [Aind2,i];
            end
        end
        num1(i,j) = length(Aind1);
        num2(i,j) = length(Aind2);
        
        %inptest
        Din = {double(D1(:,:,Aind1)), double(D2(:,:,Aind2))};
        [snr_rec(i,j,:), psnr_rec(i,j,:)] = inptest(Din,lambda_inp,S,perc_noise,maxiter);
        disp(['psnr = ', num2str(psnr_rec(i,j,:))]);
        disp(['snr = ', num2str(snr_rec(i,j,:))]);             
    end
end

save('InpTestResult.mat','mu_all','lambda_all','num1','num2','psnr_rec','snr_rec');


% %load
% load('Dict_12x12.mat');
% Din = {D};
% [a,b] = inptest(Din,.02,S,.5,300);





