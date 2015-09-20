% demo for plotting eigenvalues


tau_all = [2,4,8]; 
S = double(stdimage('lena.grey'))/255; 
S = imresize(S,.5); 
S = S(30:1:120,30:1:120); 
[Sl,Sh] = lowpass(S,4,15); 


E_all = []; 
for i = 1:length(tau_all)

    disp(num2str(i));
      
    imsz = size(Sh);
    psz = [12,12];
    stpsz = [1,1];
    %generate graph
    optl = {};
    optl.wsz = imsz;
    optl.psz = [8,8];
    optl.neig = 150;
    optl.Lformat = 'Eig';
    optl.Graph.Laplacian = 'n';
    optl.Graph.tau = tau_all(i);
    optl.Graph.Metric = 'Cosine';
    optl.Graph.GraphType = 'Full';
    optl.SaveMem = 0;
    [L,~] = laplacian_from_image(Sh,optl); 
    E_all(:,i) = sort(L{1,1}.E,'ascend');
    disp(['Graph Generated']);        
end

plot(1:150,E_all(:,1), 'r', 1:150,E_all(:,2), 'g', 1:150,E_all(:,3), 'b');
h = gca; 
xlabel('Eigenvalue No.'); 
ylabel('Eigenvalue'); 
legend('tau = 2' , 'tau = 4', 'tau = 6'); 
set(h,'fontsize',22); 


