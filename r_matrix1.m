% matrix size: 16000x16000, nnz: 1.28e+08
r_matrix_fid = fopen('r_matrix1.m.dat', 'r');
r_matrix = fscanf(r_matrix_fid,'%lf');
fclose(r_matrix_fid);
r_matrix = reshape(r_matrix, 3, round(length(r_matrix)/3));
r_matrix = spconvert(r_matrix');
