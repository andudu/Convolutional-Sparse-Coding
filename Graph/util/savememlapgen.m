function [L] = savememlapgen(s1d, opt)
tau = opt.tau;

    if strcmp(opt.GraphType,'KNearest')
        k = opt.k;
        
        nblktotal = 6*1e7;
        nwcol = size(s1d,2);
        nwrow = size(s1d,2);
        nrow_perblk = min(ceil(nblktotal/nwcol),nwrow);
        nblocks = ceil(nwrow/nrow_perblk);
        
        Iw = [];
        Jw = [];
        Rw = [];
        
        for iter = 1:nblocks
            ind = (iter-1)*nrow_perblk + 1: 1: min((iter)*nrow_perblk,nwrow);
            if strcmp(opt.Metric,'Euclidean')
                temp = sqdist(s1d,s1d(:,ind));
                temp = exp(-temp/(2*tau));
            else
                temp = cosdist(s1d,s1d(:,ind));
                temp(temp<=0) = 0; %hard thresholding negative values
                temp = exp(-abs((1./(temp+.01)-.9901).^(1.3))/tau);
            end


            for i = 1:size(temp,1)  %this maybe very slow!
                [v,is] = sort(temp(i,:),'descend');
                temp(i,:) = 0;
                temp(i,is(1:k)) = v(1:k);
            end
            temp(1:size(temp,1)+1:size(temp,1)*size(temp,1)) = 0; %zero the diagonals
            temp = sparse(temp);
            [itemp,jtemp,rtemp] = find(temp);
            itemp = itemp + (ind(1)-1);
            Iw = [Iw; itemp];
            Jw = [Jw; jtemp];
            Rw = [Rw; rtemp];          
        end 
        
        W = sparse(Iw,Jw,Rw,size(s1d,2),size(s1d,2));
        W = max(W,W');
        
        if opt.Laplacian == 'n'
           d = sqrt(sum(W,2));
           W = bsxfun(@rdivide,W,d);
           W = bsxfun(@rdivide,W,d');
           L = speye(size(W,1),size(W,2)) - W;
           
        else
           d = sum(W,2);
           L = sparse(1:1:size(W,1),1:1:size(W,1),d) - W;
        end
        
    end
    
    
    if strcmp(opt.GraphType,'WindowKNearest')  %Windowing in Large memory cases
        k = opt.k;
        
        nblktotal = 6*1e7;
        nwcol = size(s1d,2);
        nwrow = size(s1d,2);
        nrow_perblk = min(ceil(nblktotal/nwcol),nwrow);
        nblocks = ceil(nwrow/nrow_perblk);
        
        Iw = [];
        Jw = [];
        Rw = [];
        
        M = spmask(opt.coefsz(1),opt.nsz(1));
        
        for iter = 1:nblocks
            ind = (iter-1)*nrow_perblk + 1: 1: min((iter)*nrow_perblk,nwrow);
            if strcmp(opt.Metric,'Euclidean')
                temp = sqdist(s1d,s1d(:,ind));
                temp = exp(-temp/(2*tau));
            else
                temp = cosdist(s1d,s1d(:,ind));
                temp(temp<=0) = 0; %hard thresholding negative values
                temp = exp(-abs((1./(temp+.01)-.9901).^(1.3))/tau);
            end
            
            Mtemp = M(ind(1):1:ind(1)+size(temp,1)-1,:);
            temp = temp.*Mtemp;
            clear Mtemp; 
            for i = 1:size(temp,1)  %this maybe very slow!
                [v,is] = sort(temp(i,:),'descend');
                temp(i,:) = 0;
                temp(i,is(1:k)) = v(1:k);
            end
            temp(1:size(temp,1)+1:size(temp,1)*size(temp,1)) = 0; %zero the diagonals
            temp = sparse(temp);
            [itemp,jtemp,rtemp] = find(temp);
            itemp = itemp + (ind(1)-1);
            Iw = [Iw; itemp];
            Jw = [Jw; jtemp];
            Rw = [Rw; rtemp];          
        end 
        
        W = sparse(Iw,Jw,Rw,size(s1d,2),size(s1d,2));
        W = max(W,W');
        
        if opt.Laplacian == 'n'
           d = sqrt(sum(W,2));
           W = bsxfun(@rdivide,W,d);
           W = bsxfun(@rdivide,W,d');
           L = speye(size(W,1),size(W,2)) - W;
           
        else
           d = sum(W,2);
           L = sparse(1:1:size(W,1),1:1:size(W,1),d) - W;
        end
        
    end    
    
    
    
end

    


