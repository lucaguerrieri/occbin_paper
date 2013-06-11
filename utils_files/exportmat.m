function exportmat(hmat,filename)
m= size(hmat,1);
n= size(hmat,2);

hmatr= reshape(hmat,m*n,1);

fid = fopen(filename,'w');

for i = 1:m*n
    fprintf(fid,'%12.16f\n',hmatr(i));
end

fclose(fid);