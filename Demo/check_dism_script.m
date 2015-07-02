%check bd_ism

ah = rand(50,50,10,10);
d =  rand(1,1,10);
b = rand(50,50,10);

tic;
for i = 1:100

x1 = solvemdbd_ism(ah,d,b);

end
toc

tic;
for i = 1:100
x2 = solvemdbd_cg(ah,d,b,10e-5);
end
toc;

max(abs(vec(x1-x2)))
