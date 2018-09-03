function B = CalcA( f , k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interflag = 0;%0


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 512;%128;%256;%512
% % Im = double(imresize(rgb2gray(imread('Artery_Phantom.png')),[n n]))/255;
% Im = double(rgb2gray(imread('SynteticArtery2.png')))/255; Im = Im(:,337:end);Im = imresize(Im,[n n]);
% f=Im;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


interX=0;%
if interflag
    interX=3;
    f = imresize(f,interX);
end
%K=0.005;
%g = @(x) 1./(1+(x/K).^2);


% hSmooth = fspecial('gaussian', 15,2);
% hSmooth = fspecial('gaussian', 3, 0.4);% numerical exp
hSmooth = fspecial('gaussian', 3,1.5);%Artery 1%%%controls the smoothing of theimage
% hSmooth = fspecial('gaussian', 20,10);%Artery 2+noise
f = imfilter(f, hSmooth);
[fx, fy] = grad(f);
sfx = imfilter(fx, hSmooth);
sfy = imfilter(fy, hSmooth);
% hSmooth = fspecial('gaussian', 10 ,1.5);
% hSmooth = fspecial('gaussian', 3 ,1);% numerical exp
hSmooth = fspecial('gaussian', 3 ,1.5);%Artery 1 %%%controls the smoothing of the structure tensor
% hSmooth = fspecial('gaussian', 5 ,1.5);%Artery 2+noise

fx2 = imfilter(sfx.*sfx, hSmooth);
fxfy = imfilter(sfx.*sfy, hSmooth);
fy2 = imfilter(sfy.*sfy, hSmooth);


% BW = imbinarize(f);
% Bboundaries = bwboundaries(BW);
% imageE = zeros(size(f));
% idx = Bboundaries{1}(:,1);idy = Bboundaries{1}(:,2);
% for i=1:size(Bboundaries{1},1)
%     imageE(idx(i),idy(i))=1;
% end
% imageE = imageE(:);

[szy, szx]=size(f);
ST = [fx2(:) fxfy(:) fxfy(:) fy2(:)]';

% Create Stats for choosing K 
for n = 1 : size(ST,2)
    structureTensor = reshape(ST(:,n),2,2);
    [~, D] = eig(structureTensor);
    maxEigs(n) = max(D(:));
end
maxImEigs = reshape(maxEigs,szy,[]);
K = k * mean(maxEigs);
% Dref = 0.25;%0.25for C Shape
% Dref = 0.15;%0.15for Ellipse
% Dref = 0.05;%0.015for Artery1
% Dref = 0.0001;%0.0001for Artery1
% Dref = 0.001;%0.001for Artery2
Dref = 0.0005;%0.001for Artery2
% Create B
maxImD = zeros(1,size(ST,2));
NewST = zeros(size(ST));
for n = 1 : size(ST,2)
    %     NewST(:,n) = TransformST(ST(:,n),K);
    structureTensor = reshape(ST(:,n),2,2);
    [V, D] = eig(structureTensor);
    if (D(1,1) > D(2,2))
        %D(1,1) = g(absf);
        D(1,1) = WeickertG(D(1,1),K).^(0.5);
        %             D(1,1) = k;
        
        %             D(1,1) = k;
        
%         if sum(D(:))~=0, D(1,1) = k;           
%         if sum(D(:))>Dref, D(1,1) = k;
%         if imageE(n)==1, D(1,1) = k;
%         else, D(1,1) = 1; end
        maxImD(n) = D(1,1) ;
%         maxImD(n) = sum(D(:));
        D(2,2) = 1;
    else
        
        %D(2,2) = g(absf);
        D(2,2) = WeickertG(D(2,2),K).^(0.5);
        
        %             D(2,2) = k;
        
%         if sum(D(:))~=0, D(2,2) = k;
%         if sum(D(:))>Dref, D(2,2) = k;
%         if imageE(n)==1, D(2,2) = k;    
%         else, D(2,2) = 1; end
        maxImD(n) = D(2,2) ;
%         maxImD(n) = sum(D(:));
        D(1,1) = 1;
    end
    
    ST_temp = V*D/V;
    NewST(:,n) = ST_temp(:);
end

maxImD = reshape(maxImD,szy,[]);
B{1,1} = reshape(NewST(1,:),szy,[]);
B{1,2} = reshape(NewST(2,:),szy,[]);
B{2,1} = reshape(NewST(3,:),szy,[]);
B{2,2} = reshape(NewST(4,:),szy,[]);
% 

if interflag
    B{1,1}= imresize( B{1,1},1/interX);
    B{1,2}= imresize( B{1,2},1/interX);
    B{2,1}= imresize( B{2,1},1/interX);
    B{2,2}= imresize( B{2,2},1/interX);
    maxImEigs= imresize( maxImEigs,1/interX);
    maxImD= imresize( maxImD,1/interX);
    f= imresize( f,1/interX);
end

figure(201); imshow(maxImEigs,[]);
figure(202); imshow(maxImD,[]);
drawnow;pause(0.1);
end


