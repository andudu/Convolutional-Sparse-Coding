%script for testing the mini lasso solver


%random test data
n = 400; neig = 400;
temp = randn(n,n);
temp = temp'*temp;
[phi,E] = eigs(temp,neig,'sr');
E = randn(size(E,1),1);
Z = rand(size(phi,1),60);
L = phi*diag(E)*phi';

opt.Y_bar = [];
opt.V = [];
opt.MaxMainIter = 50;
opt.verbose = 1;
lambda = .1;
mu = .1;
rho = 1;

tic;
[Y_bar,V,o1] = eiglasso( phi, E, Z,lambda,mu, rho,opt);
toc;

opt1.Y = [];
opt1.eta = 1.2;
opt1.MaxMainIter = 50;
opt1.el = .5;
opt1.verbose = 1;
opt1.L1Weight = 1;
lambda = .1;

tic;
[Y, o2] = lasso_fista(L,Z,lambda,mu,rho,opt1);
toc;

norm(vec(Y-Y_bar))/norm(vec(Y))

plot(log(o1.Jfn-min(o1.Jfn)),'r');
hold on;
plot(log(o2.Jfn-min(o2.Jfn)),'b');
hold off;

