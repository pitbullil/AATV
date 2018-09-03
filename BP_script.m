clear all
close all
%load('AMAT.mat')

Result_dir= 'C:\Users\User\Documents\MATLAB\Results\'
%load sigMatTVL1
Phantom='Phantom2'
r_sensor = 39e-3;
% distance of the transducer from the centre of rotation [m]

n = 256;
fs=4e7;
c = [1478:0.1:1482];  % speed of sound [m/s]

image_width = 30e-3;                                % width of the image
dx = image_width/n;                                 % increment in x
time_res =1
%fac = fs./c;
limits(1) = r_sensor-(image_width)*sqrt(2)/2;       % limits for the signal
limits(2) = r_sensor+(image_width)*sqrt(2)/2;       % limits for the signal
%pos_start = max(1,int32((limits(1))*fac));
%pos_end = min(len,int32((limits(2))*fac));

%dt = dx/(time_res*c);                              % increment in t employed to make the model-based reconstruction
resfile = fullfile('sigMatHair.mat');
load(resfile)
len = size(sigMat,1);                               % number of samples in each projection
n_proj = size(sigMat,2);                        % number of projections
%ts = 1/fs:1/fs:len/fs;                              % sampling instants
%t = ts(pos_start):dt:ts(pos_end);
% sizeT = length(t);
% for i=1:n_proj
%     sigMat2(:,i) = interp1(ts,sigMat(:,i),t);
% end
% PPa1 = reshape(sigMat2,sizeT*n_proj,1);




level=4;
wname='db4';

filter_f = [0.05e6,7e6] ;                           % bandpass filter, lower and upper limits [Hz]
%% 270 degree top
ang_ini = -0.7853981633;
ang_end = 3.926990817079;
ang_step = 0.0184799567858;



time_res = 2;                                       % time resolution for model-based
c=1444
zz=205
for i = 1
    
sigMat1 = [zeros(zz,256) ;sigMat(200:end,:)];

    [Recon_BP]=Backprojection(sigMat1,n,c,image_width,r_sensor,fs,filter_f,ang_ini,ang_end,ang_step,time_res);
figure;
zb=image_width*100/2;

bp=imagesc([-zb zb],[-zb zb],-Recon_BP);
%axis off;
xlabel('x(cm)','fontsize',12);ylabel('y(cm)','fontsize',12);colorbar;
colormap(gray);
h = z(k)
if h-floor(h) == 0
    h = uint8(floor(h))
end
CTR = abs(max(Recon_BP(:))- min(Recon_BP(:)));
title(['BackProjection c=' num2str(c(i)) ' zz =' num2str(zz(i)) 'RNG = ' num2str(CTR)])
end
% TVWeight = 1e4; 	% Weight for TV penalty
% L1Weight = 1e5;	% Weight for Transform L1 penalty
% Itnlim = 50;		% Number of iterations
% 
% %generate transform operator
% WT = Wavelet('Daubechies',4,4);	% Wavelet
% TV = @(x)D(x);
% GTV = @(x)adjD(x);
% 
% % initialize Parameters for reconstruction
% param = init;
% param.N=n;
% param.MM = A_mat;
% param.L1 = WT;
% param.TV = TV;
% param.GTV = GTV;
% param.data = PPa1;
% param.TVWeight =TVWeight;     % TV penalty 
% param.L1Weight = L1Weight;  % L1 wavelet penalty
% param.Itnlim = Itnlim;
% param.noneg=0;
% 
% 
% res=zeros(n,n);
% %res=WT*res;
% tic;
% res = TVL1(res,param);
% Recon1_tvl1 = WT'*res;
%  tvl1_time =toc;
% Recon1_tvl1=Recon1_tvl1/max(max(Recon1_tvl1));
%   save .\wavelet_line\point_tvl1 Recon1_tvl1
% 
% 
% figure;
% point_tvl1=imagesc([-zb zb],[-zb zb],Recon1_tvl1);
% title('TV-L1')
% xlabel('x(cm)','fontsize',12);ylabel('y(cm)','fontsize',12);colorbar;
% saveas(point_tvl1,'.\wavelet_line\point_tvl1_color.jpg') 
% saveas(point_tvl1,'.\wavelet_line\point_tvl1.fig') 
% colormap(gray);
% saveas(point_tvl1,'.\wavelet_line\point_tvl1_gray.jpg') 
%  Recon_point= lsqr_b(A_mat, sparse(PPa1), 20);
% point_lsqr_time=toc(tstart)
% 
% save .\wavelet_line\point_lsqr Recon_point
% Recon1_lsqr=reshape(Recon_point(:,end),n,n);
% Recon1_lsqr=Recon1_lsqr/max(max(Recon1_lsqr));
% 
% 
% 
%  figure;
%  point_lsqr=imagesc([-zb zb],[-zb zb],Recon1_lsqr);
%  title('LSQR')
% 
%  %axis off;
% xlabel('x(cm)','fontsize',12);ylabel('y(cm)','fontsize',12);colorbar;
%  saveas(point_lsqr,'.\wavelet_line\point_lsqr_color.jpg') 
%  saveas(point_lsqr,'.\wavelet_line\point_lsqr.fig')
%  colormap(gray);
%  saveas(point_lsqr,'.\wavelet_line\point_lsqr_gray.jpg') 
% 
