function A_mat_p = Calculate_MatrixMB_Luis_Proj(x_pt,y_pt,R_pt,theta,image_width,n,n_rows,proj);

nn = n*n; % number of columns of the matrix
lt = size(x_pt,2); % length of the time vector
n_angles = size(x_pt,1); % number of points of the curve
Dxy = image_width/(n-1); % sampling distance in x and y

x_pt_unrot = x_pt*cos(theta)-y_pt*sin(theta); % horizontal position of the points of the curve in the original grid (not rotated)
y_pt_unrot = x_pt*sin(theta)+y_pt*cos(theta); % vertical position of the points of the curve in the original grid (not rotated)

x_aux_1 = [x_pt_unrot; zeros(1,lt)];
y_aux_1 = [y_pt_unrot; zeros(1,lt)];
x_aux_2 = [zeros(1,lt); x_pt_unrot];
y_aux_2 = [zeros(1,lt); y_pt_unrot];
dist_aux = sqrt((x_aux_1-x_aux_2).^2+(y_aux_1-y_aux_2).^2);
clear x_aux_1 x_aux_2 y_aux_1 y_aux_2
d_pt = dist_aux(2:n_angles,:); % length of the segments of the curve
clear dist_aux

vec_int = (1/2)*([d_pt ; zeros(1,lt)]+[zeros(1,lt) ; d_pt])./R_pt; % vector for calculating the integral
clear d_pt

x_pt_pos_aux = (x_pt_unrot+(image_width/2))/Dxy+1; % horizontal position of the points of the curve in normalized coordinates
clear x_pt_unrot
y_pt_pos_aux = (y_pt_unrot+(image_width/2))/Dxy+1; % vertical position of the points of the curve in normalized coordinates
clear y_pt_unrot

x_pt_pos_bef = floor(x_pt_pos_aux); % horizontal position of the point of the grid at the left of the point (normalized coordinates)
x_pt_pos_aft = floor(x_pt_pos_aux+1); % horizontal position of the point of the grid at the right of the point (normalized coordinates)
x_pt_dif_bef = x_pt_pos_aux-x_pt_pos_bef;
clear x_pt_pos_aux
y_pt_pos_bef = floor(y_pt_pos_aux); % vertical position of the point of the grid below of the point (normalized coordinates)
y_pt_pos_aft = floor(y_pt_pos_aux+1); % vertical position of the point of the grid above of the point (normalized coordinates)
y_pt_dif_bef = y_pt_pos_aux-y_pt_pos_bef;
clear y_pt_pos_aux

pos_triang_1x = x_pt_pos_bef; pos_triang_1y = y_pt_pos_bef; % position of the first point of the triangle
pos_triang_2x = x_pt_pos_aft; pos_triang_2y = y_pt_pos_aft; % position of the second point of the triangle
pos_triang_3ax = x_pt_pos_aft; pos_triang_3ay = y_pt_pos_bef; % position of the third point of the triangle (first case)
pos_triang_3bx = x_pt_pos_bef; pos_triang_3by = y_pt_pos_aft; % position of the third point of the triangle (second case)
pos_triang_3x = (x_pt_dif_bef>=y_pt_dif_bef).*pos_triang_3ax...
    +(x_pt_dif_bef<y_pt_dif_bef).*pos_triang_3bx; % horizontal position of the third point of the triangle
pos_triang_3y = (x_pt_dif_bef>=y_pt_dif_bef).*pos_triang_3ay...
    +(x_pt_dif_bef<y_pt_dif_bef).*pos_triang_3by; % vertical position of the third point of the triangle
clear pos_triang_3ax pos_triang_3bx pos_triang_3ay pos_triang_3by
in_pos_triang_1 = (pos_triang_1x>0)&(pos_triang_1x<=n)&(pos_triang_1y>0)&(pos_triang_1y<=n); % boolean determining if the first point of the triangle is inside the grid
in_pos_triang_2 = (pos_triang_2x>0)&(pos_triang_2x<=n)&(pos_triang_2y>0)&(pos_triang_2y<=n); % boolean determining if the second point of the triangle is inside the grid
in_pos_triang_3 = (pos_triang_3x>0)&(pos_triang_3x<=n)&(pos_triang_3y>0)&(pos_triang_3y<=n); % boolean determining if the third point of the triangle is inside the grid
Pos_triang_1_t = n*(pos_triang_1x-1)+pos_triang_1y; % one dimensional position of the first points of the triangles in the grid
Pos_triang_2_t = n*(pos_triang_2x-1)+pos_triang_2y; % one dimensional position of the second points of the triangles in the grid
Pos_triang_3_t = n*(pos_triang_3x-1)+pos_triang_3y; % one dimensional position of the third points of the triangles in the grid
clear pos_triang_1x pos_triang_1y pos_triang_2x pos_triang_2y pos_triang_3x pos_triang_3y
Pos_triang_1_t_vec = reshape(Pos_triang_1_t,1,n_angles*lt); % Pos_triang_1_t in vector form
Pos_triang_2_t_vec = reshape(Pos_triang_2_t,1,n_angles*lt); % Pos_triang_2_t in vector form
Pos_triang_3_t_vec = reshape(Pos_triang_3_t,1,n_angles*lt); % Pos_triang_3_t in vector form
clear Pos_triang_1_t Pos_triang_2_t Pos_triang_3_t

weight_triang_aux_1a = 1-x_pt_dif_bef; % weight of the first point of the triangle in the first case
weight_triang_aux_1b = 1-y_pt_dif_bef; % weight of the first point of the triangle in the second case
weight_triang_1 = ((x_pt_dif_bef>=y_pt_dif_bef).*weight_triang_aux_1a...
    +(x_pt_dif_bef<y_pt_dif_bef).*weight_triang_aux_1b).*vec_int; % weight of the first point of the triangle
clear weight_triang_aux_1a weight_triang_aux_1b
weight_triang_aux_2a = y_pt_dif_bef; % weight of the second point of the triangle in the first case
weight_triang_aux_2b = x_pt_dif_bef; % weight of the second point of the triangle in the second case
weight_triang_2 = ((x_pt_dif_bef>=y_pt_dif_bef).*weight_triang_aux_2a...
    +(x_pt_dif_bef<y_pt_dif_bef).*weight_triang_aux_2b).*vec_int; % weight of the second point of the triangle
clear weight_triang_aux_2a weight_triang_aux_2b
weight_triang_aux_3a = x_pt_dif_bef-y_pt_dif_bef; % weight of the third point of the triangle in the first case
weight_triang_aux_3b = y_pt_dif_bef-x_pt_dif_bef; % weight of the third point of the triangle in the second case
weight_triang_3 = ((x_pt_dif_bef>=y_pt_dif_bef).*weight_triang_aux_3a...
    +(x_pt_dif_bef<y_pt_dif_bef).*weight_triang_aux_3b).*vec_int; % weight of the third point of the triangle
clear weight_triang_aux_3a weight_triang_aux_3b
clear x_pt_pos_bef x_pt_pos_aft x_pt_dif_bef y_pt_pos_bef y_pt_pos_aft y_pt_dif_bef
weight_triang_1_t_vec = reshape(weight_triang_1,1,n_angles*lt); % weight_triang_1 in vector form
weight_triang_2_t_vec = reshape(weight_triang_2,1,n_angles*lt); % weight_triang_2 in vector form
weight_triang_3_t_vec = reshape(weight_triang_3,1,n_angles*lt); % weight_triang_3 in vector form
clear weight_triang_1 weight_triang_2 weight_triang_3

Row_Matrix = (linspace(1,lt,lt).'*ones(1,n_angles)).'; % rows of the sparse matrix
Row_Matrix_vec = reshape(Row_Matrix,1,n_angles*lt); % rows of the sparse matrix in vector form
clear Row_Matrix

A_mat_p = sparse([Row_Matrix_vec(in_pos_triang_1)+((proj-1)*lt) Row_Matrix_vec(in_pos_triang_2)+((proj-1)*lt) Row_Matrix_vec(in_pos_triang_3)+((proj-1)*lt)],...
    [Pos_triang_1_t_vec(in_pos_triang_1) Pos_triang_2_t_vec(in_pos_triang_2) Pos_triang_3_t_vec(in_pos_triang_3)],...
    [weight_triang_1_t_vec(in_pos_triang_1) weight_triang_2_t_vec(in_pos_triang_2) weight_triang_3_t_vec(in_pos_triang_3)],n_rows,nn);




