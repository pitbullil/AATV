clear all
close all



load sigMatTVL1 
sigMat;

n = 512;%512
xsize=[n,n];
Im = double(imresize(rgb2gray(imread('vascular2.jpg')),[n n]))/255;

c = 1527;  % speed of sound [m/s]
image_width = 25e-3;                                % width of the image
dx=image_width/(n-1);
zb=image_width*100/2;

x = (-1)*(n/2-0.5)*dx:dx:(n/2-0.5)*dx;                       
y = (-1)*(n/2-0.5)*dx:dx:(n/2-0.5)*dx; 
[X,Y] = meshgrid(x,y);

len = size(sigMat,1);                               % number of samples in each projection
n_proj = size(sigMat,2);                        % number of projections
r_sensor = 40.5e-3;                                 % distance of the transducer from the centre of rotation [m]
fs=4e7;

level=4;
wname='db4';

filter_f = [0.05e6,7e6] ;                           % bandpass filter, lower and upper limits [Hz]
ts = 1/fs:1/fs:len/fs;                              % sampling instants
%% 270 degree top
ang_ini = -0.7853981633;
ang_end = 3.926990817079;
ang_step = 4*0.0184799567858;
angle_sensor = [ang_ini:ang_step:ang_end];


n_proj=size(angle_sensor,2);
n_angles = 2*n;                                     % number of points for discretizing the curve
limits(1) = r_sensor-(image_width)*sqrt(2)/2;       % limits for the signal
limits(2) = r_sensor+(image_width)*sqrt(2)/2;       % limits for the signal
time_res = 2;                                       % time resolution for model-based
dx = image_width/n;                                 % increment in x
dt = dx/(time_res*c);                              % increment in t employed to make the model-based reconstruction


%% filter parameters
if( filter_f(2) ~= 0 )
    f_LPF = filter_f(2);
    [b_LPF,a_LPF] = cheby1( 8, .01, 2 * f_LPF/fs * .92 );
end
if( filter_f(1) ~= 0 )
    f_HPF = filter_f(1);
    if( 2 * f_HPF/fs < 0.016 )
        [b_HPF,a_HPF] = cheby1( 2, .01, 2 * f_HPF/fs * 3.3, 'high' );
    else
        [b_HPF,a_HPF] = cheby1( 4, .01, 2 * f_HPF/fs * 1.46, 'high' );
    end
end

% filtering
if( filter_f(2) ~= 0 )
    sigMat = filtfilt(b_LPF, a_LPF, sigMat);
end;
if( filter_f(1) ~= 0 )
    sigMat = filtfilt(b_HPF, a_HPF, sigMat);
end


%% back projection
image_select='full';
xc=0;
yc=0;

tstart=tic;
[Recon_BP,X,Y] = backproject_luis(sigMat,n,r_sensor,angle_sensor,c,image_select,ts,fs,image_width,xc,yc);
save ([pwd '\wavelet_line\bp'], 'Recon_BP')
bp_time=toc(tstart)
Recon_BP=-1*Recon_BP/max(max(Recon_BP));
figure;
bp=imagesc([-zb zb],[-zb zb],-Recon_BP);
%axis off;
title('BP');
xlabel('x(cm)','fontsize',12);ylabel('y(cm)','fontsize',12);colorbar;
saveas(bp,[pwd '\wavelet_line\bp_color.jpg']) 
saveas(bp,[pwd '\wavelet_line\bp.fig'])
colormap(gray);
saveas(bp,[pwd '\wavelet_line\bp_gray.jpg'])


%% Cutting and downsampling

fac = fs/c;
pos_start = max(1,int32((limits(1))*fac));
pos_end = min(len,int32((limits(2))*fac));
t = ts(pos_start):dt:ts(pos_end);
sizeT = length(t);
for i=1:n_proj
    sigMat2(:,i) = interp1(ts,sigMat(:,i),t);
end
PPa1 = reshape(sigMat2,sizeT*n_proj,1);

%% Build model matrix
PPa1_old=PPa1;
if exist(['Amat' num2str(n) '.mat'])
    load (['Amat' num2str(n) '.mat'])
else
    A_mat = Calculate_MatrixMB_Luis(c,n,image_width,t,r_sensor,angle_sensor,n_angles);
end
PPa1 = A_mat*Im(:);
PPa1_orig = PPa1;
RandomNoise = rand(size(PPa1));
PPa1_deformed = PPa1 + .8*max(PPa1)*(RandomNoise - mean(RandomNoise));

% figure; plot(PPa1); hold on;  plot(PPa1, 'R')
%% point LSQR

tstart=tic;
Recon_point= lsqr_b(A_mat, sparse(PPa1), 20);
point_lsqr_time=toc(tstart)

save ([pwd '\wavelet_line\point_lsqr'], 'Recon_point')
Recon1_lsqr=reshape(Recon_point(:,end),n,n);
Recon1_lsqr=Recon1_lsqr/max(max(Recon1_lsqr));



 figure;
 point_lsqr=imagesc([-zb zb],[-zb zb],Recon1_lsqr);
title('point lsqr');
 %axis off;
xlabel('x(cm)','fontsize',12);ylabel('y(cm)','fontsize',12);colorbar;
 saveas(point_lsqr,[pwd '\wavelet_line\point_lsqr_color.jpg']) 
 saveas(point_lsqr,[pwd '\wavelet_line\point_lsqr.fig'])
 colormap(gray);
 saveas(point_lsqr,[pwd '\wavelet_line\point_lsqr_gray.jpg']) 
 




%% point TIKHONOV

tstart=tic;
corner_point=1e6;
lambda1=corner_point;
Recon_point_tikhonov= Tikhonov_lap(A_mat,sparse(PPa1),n,lambda1,20);
point_tikhonov_time=toc(tstart)
save ([pwd '\wavelet_line\point_tikhonov'], 'Recon_point_tikhonov')
Recon1_tikhonov=reshape(Recon_point_tikhonov(:,end),n,n);
Recon1_tikhonov=Recon1_tikhonov/max(max(Recon1_tikhonov));

 figure;
 point_tikhonov=imagesc([-zb zb],[-zb zb],Recon1_tikhonov);
 title('point tikhonov');
 %axis off;
xlabel('x(cm)','fontsize',12);ylabel('y(cm)','fontsize',12);colorbar;
 saveas(point_tikhonov,[pwd '\wavelet_line\point_tikhonov_color.jpg']) 
 saveas(point_tikhonov,[pwd '\wavelet_line\point_tikhonov.fig'])
 colormap(gray);
 saveas(point_tikhonov,[pwd '\wavelet_line\point_tikhonov_gray.jpg'])

 %% TVL1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%TVL1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TVWeight = 5e10; 	% Weight for TV penalty 5e4
L1Weight = 1e5;	% Weight for Transform L1 penalty 1e5
Itnlim = 50;		% Number of iterations 50

%generate transform operator
WT = Wavelet('Daubechies',4,4);	% Wavelet
TV = @(x)D(x);
GTV = @(x)adjD(x);

% initialize Parameters for reconstruction
param = init;
param.N=n;
param.MM = A_mat;
param.L1 = WT;
param.TV = TV;
param.GTV = GTV;
param.data = PPa1;
param.TVWeight =TVWeight;     % TV penalty 
param.L1Weight = L1Weight;  % L1 wavelet penalty
param.Itnlim = Itnlim;
param.noneg=0;


res=zeros(n,n);
%res=WT*res;
tic;
res = TVL1(res,param);
Recon1_tvl1 = WT'*res;
toc;

Recon1_tvl1=Recon1_tvl1/max(max(Recon1_tvl1));
  save ([pwd '\wavelet_line\point_tvl1'], 'Recon1_tvl1')


figure;
point_tvl1=imagesc([-zb zb],[-zb zb],Recon1_tvl1);
title('point tvl1');
xlabel('x(cm)','fontsize',12);ylabel('y(cm)','fontsize',12);colorbar;
saveas(point_tvl1,[pwd '\wavelet_line\point_tvl1_color.jpg']) 
saveas(point_tvl1,[pwd '\wavelet_line\point_tvl1.fig']) 
colormap(gray);
saveas(point_tvl1,[pwd '\wavelet_line\point_tvl1_gray.jpg']) 

 %% AATVL1 _ take2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%ATVL1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ATVWeight = 5e10; 	% Weight for TV penalty 5e4
L1Weight = 1e5;	% Weight for Transform L1 penalty 1e5
Itnlim = 50;		% Number of iterations 50

%generate transform operator
WT = Wavelet('Daubechies',4,4);	% Wavelet
ATV = @(x,A)DA(x,A);
GATV = @(x,A)adjDA(x,A);

% initialize Parameters for reconstruction
param = init;
param.N=n;
param.MM = A_mat;
param.L1 = WT;
param.ATV = ATV;
param.GATV = GATV;
param.data = PPa1;
param.ATVWeight =ATVWeight;     % TV penalty 
param.L1Weight = L1Weight;  % L1 wavelet penalty
param.Itnlim = Itnlim;
param.noneg=0;
param.k=0.7;%0.05
res=zeros(n,n);
param.B = CalcB( param.L1'*res , param.k);

res=zeros(n,n);
%res=WT*res;
tic;
res = ATVL1(res,param);
Recon1_Atvl1 = WT'*res;
toc;
% Error calcuation:


Recon1_Atvl1=Recon1_Atvl1/max(max(Recon1_Atvl1));
  save ([pwd '\wavelet_line\point_Atvl1'], 'Recon1_Atvl1')

MAD_AATV = sum(abs(Recon1_Atvl1(:) - Im(:)))./(2*length(Im))
MAD_TVL1 = sum(abs(Recon1_tvl1(:) - Im(:)))./(2*length(Im))
MAD_PointDet = sum(abs(Recon1_tikhonov(:) - Im(:)))./(2*length(Im))
  
% Calculate PSNR value: TODO: added by Gilad
PSNR_AATV = PSNR(Recon1_Atvl1, Im)
PSNR_TVL1  = PSNR(Recon1_tvl1, Im)
PSNR_PointDet = PSNR(Recon1_tikhonov, Im)



  
figure;
point_Atvl1=imagesc([-zb zb],[-zb zb],Recon1_Atvl1);
title('point Atvl1');
xlabel('x(cm)','fontsize',12);ylabel('y(cm)','fontsize',12);colorbar;
saveas(point_Atvl1,[pwd '\wavelet_line\point_Atvl1_color.jpg']) 
saveas(point_Atvl1,[pwd '\wavelet_line\point_Atvl1.fig']) 
colormap(gray);
saveas(point_Atvl1,[pwd '\wavelet_line\point_Atvl1_gray.jpg']) 
%%
% save Data512.mat -v7.3