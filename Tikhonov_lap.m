function recon= Tikhonov_lap(A_mat,b_vec,n,lambda,lsqr_iter);

nn = n*n;

FILT = lambda*calculate_matrix_highpass_4(n);

L = sparse(FILT);

cell_mat{1,1} = A_mat;
cell_mat{2,1} = L;
A_mat = cell2mat(cell_mat);
b_vec = [b_vec; zeros(nn,1)];
 recon=lsqr_b(A_mat,sparse(b_vec),lsqr_iter);

end


