
%generate

D0 = zeros(12,12,30, 'single');
D0(4:9,4:9,:) = single(randn(6,6,30));

% num_all = [5,10,20];
% lambda_all = [.05,.15,.3];
% 
% for i = 2:3
%     for j = 1:3
%         disp([num2str(i), ',', num2str(j)]);
%         n = num_all(j);
%         D_init = D0(:,:,1:n);
%         lambda = lambda_all(i);
%         mu = lambda/2;
%         fname = (['Experiment1/DictComp',num2str(i),num2str(j),'.mat']);
%         [D1,X1,D2,X2] = comptest(D_init,lambda,mu);
%         save(fname,'D1','X1','D2','X2','lambda','mu');
%     end
% end


%num_all = [6,10,16];
lambda_all = [.2,.3,.4];
mu_all = [.2,.3,.4];

for i = 1:3
    for j = 1:3
        disp([num2str(i), ',', num2str(j)]);
        n = 8;
        D_init = D0(:,:,1:n);
        lambda = lambda_all(j);
        mu = mu_all(i);
        fname = (['Experiment3/DictComp',num2str(i),num2str(j),'.mat']);
        [D1,X1,D2,X2] = comptest(D_init,lambda,mu);
        save(fname,'D1','X1','D2','X2','lambda','mu');
    end
end
