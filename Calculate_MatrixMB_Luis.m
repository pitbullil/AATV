function A_mat = Calculate_MatrixMB_Luis(c0,n,image_width,t,r_sensor,angle_sensor,n_angles);

nn = n*n; % number of columns of the matrix
n_rows = length(t)*length(angle_sensor); % number of rows of the matrix
Dxy = image_width/(n-1);
dt = 1e-15; % diferential of time to perform derivation
tpdt = t+dt; % time instants for t+dt
t = t-dt;

A_mat = sparse([],[],[],n_rows,nn);

angle_max = asin(((image_width+2*Dxy)*sqrt(2))/(2*r_sensor));
dist_sensor = c0*t;
dist_pdd_sensor = c0*tpdt;
angles = linspace(-angle_max,angle_max,n_angles).'*ones(1,length(t));

for i=1:length(angle_sensor)
    
    disp(['Projection Number: ' num2str(i)]);
    theta = angle_sensor(i);
    
    x_pt = r_sensor-(ones(n_angles,1)*dist_sensor).*cos(angles);
    y_pt = (ones(n_angles,1)*dist_sensor).*sin(angles);
    R_pt = ones(n_angles,1)*dist_sensor;
    xpdx_pt = r_sensor-(ones(n_angles,1)*dist_pdd_sensor).*cos(angles);
    ypdy_pt = (ones(n_angles,1)*dist_pdd_sensor).*sin(angles);
    RpdR_pt = ones(n_angles,1)*dist_pdd_sensor;
    A_mat = A_mat + (1/(2*dt))*...
        (-Calculate_MatrixMB_Luis_Proj(x_pt,y_pt,R_pt,theta,image_width,n,n_rows,i)+...
        Calculate_MatrixMB_Luis_Proj(xpdx_pt,ypdy_pt,RpdR_pt,theta,image_width,n,n_rows,i));
    clear x_pt y_pt R_pt xpdx_pt ypdy_pt RpdR_pt
    
end



