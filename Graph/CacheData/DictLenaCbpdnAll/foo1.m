% extract all regular dictionaries 

o1.grey = 1;
o2.sparse = 1;

for i = 1:5
    load(['DictComp5',num2str(i),'6']);
    square_plot(D2,o1);
    square_plot(X2,o2);
end