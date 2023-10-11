% matrix size: 2000x2000, nnz: 2.0512e+06
A_fid = fopen('A1.m.dat', 'r');
A = fscanf(A_fid,'%lf');
fclose(A_fid);
A = reshape(A, 3, round(length(A)/3));
A = spconvert(A');
