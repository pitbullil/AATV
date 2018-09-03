function FILT = calculate_matrix_highpass_4(n);

pos = linspace(1,n*n,n*n).';
values_diag = (4/5)*ones(n*n,1);
FILT = sparse(pos,pos,values_diag,n*n,n*n);

pos_left = pos-n;
pos_right = pos+n;
pos_down = pos-1;
pos_up = pos+1;
cond_left_in = pos_left>0;
cond_right_in = pos_right<=n*n;
cond_down_in = floor((pos_down-1)/n)==floor((pos-1)/n);
cond_up_in = floor((pos_up-1)/n)==floor((pos-1)/n);

columns = [pos_left, pos_right, pos_down, pos_up];
rows = pos*ones(1,4);
values = (-1/5)*ones(n*n,4);
cond = [(cond_left_in), (cond_right_in), (cond_down_in), (cond_up_in) ];

FILT = FILT + sparse(rows(cond),columns(cond),values(cond),n*n,n*n);


