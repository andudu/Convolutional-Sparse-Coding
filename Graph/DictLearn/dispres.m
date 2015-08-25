function dispres(i,j)
load(['DictComp',num2str(i),num2str(j),'.mat']);
o.grey = 1;
g.sparse = 1;
square_plot(D1,o);
square_plot(X1,g);
square_plot(D2,o);
square_plot(X2,g);

end

