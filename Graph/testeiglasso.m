%script for testing the block mini lasso solver


%random test data
n = 500; neig = 50;
temp = randn(n,n);
temp = temp'*temp;
[phi,E] = eigs(temp,neig,'sr');
E = randn(size(E,1),1);

opt.Y_bar = [];
opt.V = [];
opt.Maxiter = 50;
opt.verbose = 1;
Z = rand(size(phi,1),50);
lambda = .1;
mu = .1;
rho = 1;

Y_bar = eiglasso( phi, E, Z,lambda,mu, rho,opt);

Y2 = shrink(Z,lambda/rho);

norm(vec(Y2-Y_bar))/norm(vec(Y2))


