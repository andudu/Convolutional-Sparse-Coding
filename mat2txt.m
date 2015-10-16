function mat2txt(M,fname)

fid = fopen(fname,'w');
for i = 1:size(M,1)
    for j = 1:size(M,2)
        if j == size(M,2)
            sep = '\n'; 
        else
            sep = ',';
        end
        fprintf(fid,strcat('%.4f',sep),M(i,j));
    end
end
fclose(fid);



