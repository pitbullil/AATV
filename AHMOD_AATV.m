function [u,Err,postsets]=AHMOD_AATV(A_Model,b,lambda,params,presets)
%% solves min F(Kx)+G(x)
%         x
% for:
%   G(x) = (lambda/2)*||A_Model*x-b||^2_2
%   F(Dx) = ||Dx||_1 (=TV(x))
% using the AHMOD variant described in Chambolle-Pock
%% input: A - model matrix; b - sound wave; lambda - ratio between energies;
%% output: u - processed image; Err - error in relation to the prev u.
%% example: u = AHMOD(A,b,lambda)
% %% Or,
% load('D:\Cloud\OneDrive - Technion\OptoAcousticsAATV\TVL1 for Amir\Amat128.mat')
% image = (double(rgb2gray(imread('vascular2.jpg'))))/255;
% image = imresize(image,[128,128]);
% image = image-min(image(:));
% image = image/max(image(:));
% image = image - mean(image(:));
% figure(7);imshow(image,[]);Im = image
% params.maxIter = 1000;params.print = 1; params.n = 128;A = A_mat;b = A_mat*Im(:); lambda = 8;[u,Err,postsets]=AHMOD_TV(A,b,lambda,params);

if ~exist('params') || (exist('params') && ~isfield(params,'maxIter')) params.maxIter = 1000; end
if ~exist('params') || (exist('params') && ~isfield(params,'print'))params.print = 0; end
if ~exist('params') || (exist('params') && ~isfield(params,'A')) params.A = CalcA( zeros(params.n) ,params.k);end
if ~exist('params') || (exist('params') && ~isfield(params,'CalcAInside')), params.CalcAInside = 0; end

% L_A = sqrt(norm(A,1)*norm(A,inf));
% A = A/L_A;

ProxFSz = @ (z,sigma) z./max(1,abs(z));
ProxFSq = @ (q,sigma) (q-sigma*b)/(1+sigma/lambda);
ProxFS = @(q,z,sigma) deal(ProxFSq(q,sigma), ProxFSz(z,sigma));
ProxG = @ (u,tau) u;


K_A = @ (u) A_Model*u(:);
K_nabla = @ (u) gradA(params.A, u);
KS_A = @ (q) reshape(A_Model'*q,params.n,params.n);
KS_nabla = @ (z) - divA(params.A, z);
KS = @ (q, z) KS_A(q) + KS_nabla(z);

absV = @(v) sqrt(sum(v.^2,3)); % gradient magnitude
G = @(u) (lambda/2)*norm(A_Model*u(:)-b);
F = @(v) norm(absV(v));


% L_A = 1;
L_A = sqrt(norm(A_Model,1)*norm(A_Model,inf));L_nabla = 8; L = L_A + L_nabla;
gamma = 0.7*lambda;

if ~exist('presets')
    u = zeros(params.n); % -> u = restored image
    q = zeros(size(b)); % -> q = lambda*(Ax-b)
    z = zeros(params.n,params.n,2); % -> z = xi_Omega
    tau = 0.5;%0.02article|0.05 looks good
    sigma = 1/(tau*L^2);%should satisfy sigma*tau*norm(K)^2<1
else
    u = presets.u;
    q = presets.q; 
    z = presets.z; 
    tau = presets.tau; 
    sigma = presets.sigma;     
end
    
objective = [];

for i=1:params.maxIter
    u_prev = u;
    [q, z] = ProxFS(q + sigma*K_A(u), z + sigma*K_nabla(u), sigma);
    u = ProxG(u - tau*KS(q,z),tau);
    theta = 1/sqrt(1 + 2*gamma*tau);
    tau = theta*tau;
    sigma = sigma/theta;
    if params.CalcAInside
        params.A = CalcA( u,params.k);
        K_nabla = @ (u) gradA(params.A, u);
        KS_nabla = @ (z) - divA(params.A, z);
    end
    if params.print
        % track convergence
        objectiveF(i) = F(K_nabla(u));
        objectiveG(i) = G(u);
        objective(i) = objectiveF(i) + objectiveG(i);
        if ~(i==params.maxIter)
            fprintf('iter %d obj.F(K(u)) = %.6f. G(u)= %.6f. J(u)=%.6f\n',i,F(K_nabla(u)),G(u),objective(i));
            if (mod(i,25)==0) 
                figure(989); imshow(u,[]);
                figure(990); plot(1:length(objective),objective,'k',1:length(objectiveF),objectiveF,'r',1:length(objectiveG),objectiveG,'b');
                legend('objective','F','G')
                drawnow;
            end
        end
    end
end
Err = norm(u-u_prev)/norm(u);
if params.print
    disp(['Error = ' num2str(norm(u-u_prev)/norm(u)) ]);
    figure(999); plot(1:length(objective),objective,'k',1:length(objectiveF),objectiveF,'r',1:length(objectiveG),objectiveG,'b');
    legend('objective','F','G')
end

% for speeding up sequential calling
postsets.u = u;
postsets.q = q;
postsets.z = z;
postsets.tau = tau;
postsets.sigma = sigma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions              %%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradient (forward difference)
function dP = gradA(A, P)
	fx = P(:,[2:end end])-P;
	fy = P([2:end end],:)-P;
    
    Afx = A{1,1}.*fx + A{1,2}.*fy;
    Afy = A{2,1}.*fx + A{2,2}.*fy;
    
    dP = cat(3,Afx,Afy);
%%%%%%%%%%%%%%%%%%%%%%%%%%% Divergence (backward difference)
function divP = divA(A, P)
    Px = P(:,:,1); Py = P(:,:,2);
    
    APx = A{1,1}.*Px + A{2,1}.*Py;
    APy = A{1,2}.*Px + A{2,2}.*Py;
    
	fx = APx-APx(:,[1 1:end-1]); fx(:,1) = APx(:,1); fx(:,end) = -APx(:,end-1);
	fy = APy-APy([1 1:end-1],:); fy(1,:) = APy(1,:); fy(end,:) = -APy(end-1,:);
    divP = fx+fy;