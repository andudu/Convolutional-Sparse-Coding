

x1 = bsxfun(@plus, -[3;0] , randn(2,500));
x2 = bsxfun(@plus, [3;0] , randn(2,500));
X = [x1,x2];

%plot(X(1,:),X(2,:), '.');
% opt.Metric = 'Euclidean';
% opt.GraphType = 'Full';
% opt.tau = 1;
% W = graphgen(X,opt);
% 
% L = nlap(W);
% 
% [phi,E] = eigs(L,3,'sr');
% 
% y = [phi(:,2)';phi(:,3)'];
% 
% figure;
% plot(y(1,:),y(2,:),'r.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = {};
opt.simMat = 'full';
opt.L = 'sym';
[ evec, eval, isolatedIndex, opts ] = spectralDiffusion( X',3,opt);

figure;
plot(evec(:,1),evec(:,2),'b.');