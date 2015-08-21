function [ evec, eval, isolatedIndex, opts ] = spectralDiffusion( points,numEigs,opts)
% [ evec, eval, opts ] = spectralCluster( points,numEigs,num_clust,opts)
% input parameters:
%   points      matrix of data points. Rows correspond to points, cols to
%                   dimensions.
%   simMat     distance matrix of a graph. Like adjacency matrix, but
%                   edges are some notion of distance, and unconnceted
%                   edges have weight inf.
%   numEigs     number of eigenvectors to use as coordinates in 
%                   clustering
%   num_clu
%   
%   opts        struct:     set of options. If not inputted, defaults are used
%
%       [.simMat]   full, eps nbhd, kNN, mutual kNN
%       [.weight]   >0  def = 1 scaling of the weight matrix
%       [.Lp]       [1,inf]. use Lp norm for computing distances. 
%       [.eps]
%       [.kNN]                      
%       [.L]        'unnormalized'(def)    laplace operator scaling
%                   'sym'
%                   'rw'
%        
%
%% construct weight matrix

%defaults used (calculated during compilation)

%opts.simMat = 'eps nbhd'
%opts.weight = 1;
%opts.eps = removes 80 to 85% of connections;
%opts.kNN = ceil(.1*numPoints);
%opts.L = 'unnormalized'

% check opts and set defaults. options specific to simMat or L are handled
% later
if nargin == 2
    opts = struct;
end
if nargin == 1 || numEigs == inf
    numEigs = size(points,1)-1;
end
    

if isfield(opts,'simMat')==0 || strcmp(opts.simMat,'def')
    opts.simMat = 'kNN';
elseif any(strcmp(opts.simMat,{'full';'eps nbhd';'kNN';'mutual kNN'}))==0
    error('opts.simMat = %s not recognized',opts.simMat)
end
if isfield(opts,'L') == 0 || strcmp(opts.L,'def')
    opts.L = 'unnormalized';
elseif any(strcmp(opts.L,{'unnormalized';'sym';'rw'}))==0
    error('opts.L = %s not recognized',opts.L);
end
if  isfield(opts,'Lp')==0 || strcmp(opts.Lp,'def')
    opts.Lp = 2;
elseif opts.Lp < 1
    error('opts.Lp should be in [1,inf]')
end


%construst similarity matrix simMat
%problems to fix: what if two points are equal?
%                 full case: very sensitive to opts.weight (sigma), even
%                 with data [rand(100,2)-3; randn(100,2)+2]. I don't know
%                 of a good way to pick sigma.

numPoints = size(points,1);
simMat = zeros(numPoints);
switch opts.simMat
    case 'full'
        simMat   = pdist2(points,points,'minkowski',opts.Lp);
        if isfield(opts,'weight') == 0 || strcmp(opts.weight,'def')
            %the larger between 1 and the smallest off diagonal element of
            %distance matrix
            opts.weight = min(simMat(triu(true(numPoints),1)))/.021;
        end
        simMat = exp(-(simMat/(sqrt(2)*opts.weight)).^2);
        simMat(1:numPoints+1:numPoints*numPoints) = 0;
    case 'eps nbhd'
        simMat = pdist2(points,points,'minkowski',opts.Lp);
        
        %default: find eps so that 80% to 85% of distances are ignored
        if isfield(opts,'eps')==0 || strcmp(opts.eps,'def')
            %opts.eps = .6*max(max(simMat));
            apple=simMat(triu(true(numPoints),1));
            banana=numel(apple);
            
            opts.eps = .1*max(apple);   %in general (epsup + epslow)/2
            delta = opts.eps;           %in general (epsup - epslow)/2
            iter = 1;
            f = nnz(apple(apple<opts.eps))/banana;  %f needs to be increasing in epsilon
            while (f > .2 || f < .15) && iter < 10
                delta = delta/2;
                if f >.2
                    opts.eps = opts.eps - delta;
                elseif f <.15
                    opts.eps = opts.eps + delta;
                end
                f = nnz(apple(apple<opts.eps))/banana;
                iter = iter+1;
            end
        end
        
        %default 2: taken from dbscan, may be better
        %k = .05*size(points,1);
        %[m,n]=size(points);
        %Eps=((prod(max(points)-min(points))*k*gamma(.5*n+1))/(m*sqrt(pi.^n))).^(1/n);

            simMat(simMat > opts.eps) = 0;
            simMat = simMat ~= 0;   
    case 'kNN'
        if isfield(opts,'kNN')==0 || strcmp(opts.kNN,'def')
            opts.kNN = ceil(.2*numPoints);
        end
        NNindex = kNN(points,points,opts.kNN+1);
        NNindex = NNindex(:,2:end);
        for ii = 1:numPoints
            simMat(ii,NNindex(ii,:))=1;
        end
        simMat = bsxfun(@or,simMat,simMat');
    case 'mutual kNN'
        if isfield(opts,'kNN')==0 || strcmp(opts.kNN,'def')
            opts.kNN = ceil(.1*numPoints);
        end
        NNindex = kNN(points,points,opts.kNN+1);
        NNindex = NNindex(:,2:end);
        for ii = 1:numPoints
            simMat(ii,NNindex(ii,:))=1;
        end
        simMat = bsxfun(@and,simMat,simMat');       
end

%can i replace this by only computing smallest eigs?
%construct the laplace operator and find eigs
switch opts.L
    case 'unnormalized'
        D = sum(simMat);
        isolatedIndex = find(D == 0);
        D(isolatedIndex)=[];
        simMat(:,isolatedIndex)=[];
        simMat(isolatedIndex,:)=[];
        
        L = diag(D)-simMat;
        [evec,eval] = eig(L);
    case 'sym'
        D = sum(simMat);
        isolatedIndex = find(D == 0);
        D(isolatedIndex)=[];
        simMat(:,isolatedIndex)=[];
        simMat(isolatedIndex,:)=[];
        
        sqrtDinv = diag( (1./D).^.5 );
        Lsym = eye(numPoints-length(isolatedIndex))-sqrtDinv*simMat*sqrtDinv;
        [evec,eval] = eig(Lsym);
        evec = bsxfun(@rdivide,evec, sum(evec.^2,2));
    case 'rw'
        D = sum(simMat);
        isolatedIndex = find(D == 0);
        D(isolatedIndex)=[];
        simMat(:,isolatedIndex)=[];
        simMat(isolatedIndex,:)=[];
        
        L = diag(D)-simMat;
        [evec,eval] = eig(L,diag(D));
end
eval=diag(eval);
[eval,apple] = sort(eval);
eval = eval(2:numEigs+1);
evec = evec(:,apple(2:numEigs+1));

end