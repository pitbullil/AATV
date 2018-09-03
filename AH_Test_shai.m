
%% load image
n = 256;
% Im = double(imresize(rgb2gray(imread('Artery_Phantom.png')),[n n]))/255;
% Im = double(rgb2gray(imread('SynteticArtery2.png')))/255; Im = Im(:,337:end);Im = imresize(Im,[n n]);
Im = double(rgb2gray(imread('vascular2.jpg')))/255; Im = Im(:,337:end);Im = imresize(Im,[n n]);
% Im = double(rgb2gray(imread('Moon.png')))/255; Im = imresize(Im,[n n]);

Im = Im - min(Im(:));
Im = Im/max(Im(:));
Im = Im - mean(Im(:));

figure(1); imshow(Im,[]);

%% create pressure wave


if exist(['Amat' num2str(n) '.mat'])
    load (['Amat' num2str(n) '.mat'])
else
    A_mat = Calculate_MatrixMB_Luis(c,n,image_width,t,r_sensor,angle_sensor,n_angles);
end
L_A = sqrt(norm(A_mat,1)*norm(A_mat,inf))/(8*20);%80;
A_mat = A_mat/L_A;
% A_mat = A_mat/1e5;
PPa1 = A_mat*Im(:);
PPa1_orig = PPa1;

% 1. noise in the image plain  
% RandomNoise = randn(size(Im));
% Im_deformed = Im + 1.4*max(Im(:))*(RandomNoise - mean(RandomNoise(:)));
% PPa1_deformed = A_mat*Im_deformed(:);
% PPa1 = PPa1_deformed;

% % 2. noise in the sound plain  
% RandVec = rand(size(PPa1)); RandomNoise = 0.6*max(PPa1)*(RandVec - mean(RandVec));
% PPa1_deformed = PPa1 + RandomNoise;
% PPa1 = PPa1_deformed;

% % 3. arcs in the sound plain 
% t=1:1500; 
% PPa1_projection_noise = PPa1_orig;
% PPa1_projection_noise(1:length(t)) = 10*PPa1_orig(1:length(t));
% PPa1 = PPa1_projection_noise;


% PPa1_projection_noise = PPa1_orig;
% idxS = fix(rand(1,10)*length(PPa1_projection_noise))';
% idxE = idxS+1500;
% for i=1:10
%    idxT(:,i) =  idxS(i):idxE(i);
% end
% idxT = idxT(:);
% PPa1_projection_noise(idxT) = 5*PPa1_orig(idxT);
% PPa1 = PPa1_projection_noise;

close(figure(2));figure(2);   plot(PPa1, 'r'); hold on;plot(PPa1_orig,'b'); legend('Noised','Original');title('PPa')
%% reconstruct
NoiseType = 'arcs';
NoiseType = 'rings';
NoiseType = 'noise';
%%
load (['.\Result_mat\256_' NoiseType '\PPa1.mat'],'PPa1');

lambda_vec = [0.05 0.1 0.3 0.5 1 1.5 2 2.5 3]; 

% lambda_vec = [0.0001 0.001 0.01 0.1 1 5 10 50];
for i=1:length(lambda_vec)
    params.maxIter = 5000; params.print = 1; params.n = n; A = A_mat; b = PPa1;
    lambda = lambda_vec(i);
    % tic;[u,Err,postsets]=AHMOD_TV(A,b,lambda,params);toc;
    params.k = 1; params.CalcAInside = 1;
    tic;[u,Err,postsets]=AHMOD_AATV(A,b,lambda,params);toc;
    % tic;[u,Err,postsets]=AHMOD_TV(A,b,lambda,params,postsets);toc;
    figure;imshow(u,[]);
    XTV.lambda = lambda;
    XTV.A = A; XTV.b = b;
    XTV.PPa1_orig = PPa1_orig;
    % XTV.RandomNoise = RandomNoise;
    XTV.u = u; XTV.Im = Im;
    XTV.postsets = postsets;
    Data_AATV{i}.XTV = XTV;
    % save('AppTest10_AHMOD_Vascular2_RandomArcSoundPlainNoise_10arcs_5mult__Lambda_0_100_TV_z_en.mat','XTV')
    % save(['AppTest2_AHMOD_Vascular2_ArcSound_5__Lambda_' num2str(lambda_vec(i)*10) '_TV_z_en.mat'],'XTV')
end

% lambda_vec = [0.0001 0.001 0.01 0.1 1 5 10 50];
for i=1:length(lambda_vec)
    params.maxIter = 5000; params.print = 1; params.n = n; A = A_mat; b = PPa1;
    lambda = lambda_vec(i);
    tic;[u,Err,postsets]=AHMOD_TV(A,b,lambda,params);toc;
    % tic;[u,Err,postsets]=AHMOD_TV(A,b,lambda,params,postsets);toc;
    figure;imshow(u,[]);
    XTV.lambda = lambda;
    XTV.A = A; XTV.b = b;
    XTV.PPa1_orig = PPa1_orig;
    % XTV.RandomNoise = RandomNoise;
    XTV.u = u; XTV.Im = Im;
    XTV.postsets = postsets;
    % save('AppTest10_AHMOD_Vascular2_RandomArcSoundPlainNoise_10arcs_5mult__Lambda_0_100_TV_z_en.mat','XTV')
    % save(['AppTest2_AHMOD_Vascular2_ArcSound_5__Lambda_' num2str(lambda_vec(i)*10) '_TV_z_en.mat'],'XTV')
    Data_TV{i}.XTV = XTV;
end
save(['Data_' NoiseType '3.mat'],'Data_TV','Data_AATV','PPa1','-v7.3')
% save data_noise.mat
% XTV_TV=XTV;
%%
figure;
load Arterymask256.mat;

DataX = Data_TV;
nTests = length(DataX);
for iTest = 1:nTests
    u = DataX{iTest}.XTV.u; lambda = DataX{iTest}.XTV.lambda; Im = DataX{iTest}.XTV.Im;
    MAD = CalcMAD(u,Im); CNR = calc_contrast(u,Bin);
    subplot(2,nTests+1,iTest); imagesc(u); title({[' \lambda = ' num2str(lambda)],['. MAD = ' num2str(MAD) 10 '. CNR = ' num2str(CNR)]});
end
DataX = Data_AATV;
nTests = length(DataX);
for iTest = 1:nTests
    u = DataX{iTest}.XTV.u; lambda = DataX{iTest}.XTV.lambda; Im = DataX{iTest}.XTV.Im;
    MAD = CalcMAD(u,Im); CNR = calc_contrast(u,Bin);
    subplot(2,nTests+1,iTest+nTests+1); imagesc(u); title({[' \lambda = ' num2str(lambda)],['. MAD = ' num2str(MAD) 10 '. CNR = ' num2str(CNR)]});
end
load (['.\Result_mat\256_' NoiseType '\TVL1_LSQR.mat']);
u = Recon1_lsqr;MAD = CalcMAD(u,Im); CNR = calc_contrast(u,Bin);
subplot(2,nTests+1,nTests+1); imagesc(u); title({['LSQR'],['. MAD = ' num2str(MAD) 10 '. CNR = ' num2str(CNR)]});
u = Recon1_tvl1;MAD = CalcMAD(u,Im); CNR = calc_contrast(u,Bin);
subplot(2,nTests+1,2*nTests+2); imagesc(u); title({['TVL1'],['. MAD = ' num2str(MAD) 10 '. CNR = ' num2str(CNR)]});
%%
drawnow;close(figure(100));figure(100);
set(gcf,'units','normalized','position',[0.1 0.1 0.35 0.44]);
[ha, pos] = tight_subplot(2,3,[.001 .001],[.02 .07],[.02 .02]) ;
set(gcf,'DefaultAxesFontSize',14)

rect1 = [143 1 50 50];rect2 = [200 90 50 50];rect3 = [5 200 50 50];rect4 = [100 150 50 50];
% set(gca, 'Position',pos{5} + [0 0.04 0 0]);

InsertRects = @(lw) [rectangle('Position',rect1,'EdgeColor','r','LineWidth',lw);...
                rectangle('Position',rect2,'EdgeColor','g','LineWidth',lw);...
                rectangle('Position',rect3,'EdgeColor','b','LineWidth',lw);...
                rectangle('Position',rect4,'EdgeColor','m','LineWidth',lw)];          
GetRectIm = @(Im,rect) Im(rect(2):rect(2)+rect(4)-1,rect(1):rect(1)+rect(3)-1);
AssembleRectIm =  @(Im) [GetRectIm(Im,rect1) GetRectIm(Im,rect2); GetRectIm(Im,rect3) GetRectIm(Im,rect4)];
lim =[];%[-0.6 0.8];% [-0.9 1.2];
axes(ha(1)); u = Recon1_lsqr; imshow(u,lim);InsertRects(1);title('LSQR','interpreter','latex');
axes(ha(1+3*1)); u = AssembleRectIm(Recon1_lsqr); imshow(u,lim);
axes(ha(2)); u = Recon1_tvl1; imshow(u,lim);InsertRects(1);title('TVL1','interpreter','latex');
axes(ha(2+3*1)); u = AssembleRectIm(Recon1_tvl1); imshow(u,lim);
Exp = 7;uAATV = Data_AATV{Exp}.XTV.u;lambda = Data_AATV{Exp}.XTV.lambda;
axes(ha(3)); u = uAATV; imshow(u,lim);InsertRects(1);title('AATV','interpreter','latex');
axes(ha(3+3*1)); u = AssembleRectIm(uAATV); imshow(u,lim);
AATVparams.nIter = params.maxIter;AATVparams.lambda = lambda;
save (['.\Result_mat\256_' NoiseType '\AATV.mat'],'uAATV', 'AATVparams');

u = Recon1_lsqr;MAD = CalcMAD(u,Im); CNR = calc_contrast(u,Bin);
disp(['LSQR: MAD = ' num2str(MAD)  '. CNR = ' num2str(CNR)])
u = Recon1_tvl1;MAD = CalcMAD(u,Im); CNR = calc_contrast(u,Bin);
disp(['TVL1: MAD = ' num2str(MAD)  '. CNR = ' num2str(CNR)])
u = Data_AATV{Exp}.XTV.u;MAD = CalcMAD(u,Im); CNR = calc_contrast(u,Bin);
disp([NoiseType ': MAD = ' num2str(MAD)  '. CNR = ' num2str(CNR)])

if ~(exist('AImages','dir') == 7), mkdir('AImages'); end
if ~(exist('AFigs','dir') == 7), mkdir('AFigs'); end
saveas(figure(100),[pwd '\AImages\Optoacoustic_Hans_Comparison_' NoiseType '_Noise_Example.png'])
saveas(figure(100),[pwd '\AFigs\Optoacoustic_Hans_Comparison_' NoiseType '_Noise_Example.fig']) 

%%
Data = Data_TV; zoomxy = [25 25+50 130 130+50];
X_TV.Zx = flipud(fliplr(Data{2}.XTV.postsets.z(:,:,1)));X_TV.Zy = flipud(fliplr(Data{2}.XTV.postsets.z(:,:,2)));
divZ = - div(flipud(fliplr(Data{2}.XTV.postsets.z)));
X_TVzoom.Zx = flipud(fliplr(Data{2}.XTV.postsets.z(zoomxy(1):zoomxy(2),zoomxy(3):zoomxy(4),1)));X_TVzoom.Zy = flipud(fliplr(Data{2}.XTV.postsets.z(zoomxy(1):zoomxy(2),zoomxy(3):zoomxy(4),2)));
Data = Data_AATV;
X_AATV.Zx = flipud(fliplr(Data{2}.XTV.postsets.z(:,:,1)));X_AATV.Zy = flipud(fliplr(Data{2}.XTV.postsets.z(:,:,2)));
X_AATVzoom.Zx = flipud(fliplr(Data{2}.XTV.postsets.z(zoomxy(1):zoomxy(2),zoomxy(3):zoomxy(4),1)));X_AATVzoom.Zy = flipud(fliplr(Data{2}.XTV.postsets.z(zoomxy(1):zoomxy(2),zoomxy(3):zoomxy(4),2)));
divAZ = - flipud(fliplr(divA(CalcA(Data{2}.XTV.u,1) ,Data{2}.XTV.postsets.z)));
figure(256);
subplot 221;CQuiver( X_TV, 1 ); subplot 223;CQuiver( X_TVzoom, 1 );
subplot 222;CQuiver( X_AATV, 1 ); subplot 224;CQuiver( X_AATVzoom, 1 );
figure(301);h1 = CQuiver( X_TV, 1 );
figure(302);h2 = CQuiver( X_TVzoom, 1 );
figure(303);h3 = CQuiver( X_AATV, 1 );
figure(304);h4 = CQuiver( X_AATVzoom, 1 );
saveas(h1,[pwd '\AImages\Optoacoustic_EclipsedMoon_gradientFlow_TV_Example.png'])
saveas(h2,[pwd '\AImages\Optoacoustic_EclipsedMoon_gradientFlow_TVzoom_Example.png'])
saveas(h3,[pwd '\AImages\Optoacoustic_EclipsedMoon_gradientFlow_AATV_Example.png'])
saveas(h4,[pwd '\AImages\Optoacoustic_EclipsedMoon_gradientFlow_AATVzoom_Example.png'])

saveas(h1,[pwd '\AFigs\Optoacoustic_EclipsedMoon_gradientFlow_TV_Example.fig'])
saveas(h2,[pwd '\AFigs\Optoacoustic_EclipsedMoon_gradientFlow_TVzoom_Example.fig'])
saveas(h3,[pwd '\AFigs\Optoacoustic_EclipsedMoon_gradientFlow_AATV_Example.fig'])
saveas(h4,[pwd '\AFigs\Optoacoustic_EclipsedMoon_gradientFlow_AATVzoom_Example.fig'])


%%
u = XTV.u;
MAD_AATV = sum(abs(u(:) - Im(:)))./(2*length(Im))
figure;
subplot 121; imshow(u,[]);
subplot 122; imshow(Im - u,[]);title(num2str(MAD_AATV));
%%
err1 = -postsets.q/lambda;
err2 = (0.4*max(PPa1)*(RandomNoise - mean(RandomNoise)));
close(figure(5));figure(5);   plot(err1, 'r'); hold on;plot(err2,'b'); legend('Recons Noise','Original Noise');title('Noise')
close(figure(6));figure(6);   plot(err1-err2); title('q-noise')
