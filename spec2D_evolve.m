function XTV = spec2D_evolve(b, params)
% function XTV = spec2D_evolve(b, params)
% based on a private function by Guy Gilboa
% Computes TV/AATV spectrum of an image, returns S, time interval T, Phi(T)
% and residual image f_r (to be added in the reconstruction).
% Example: XTV = spec2D_evolve(b) - using default params: TV regularizer based reconstruction.
% Based on: [1] G. Gilboa, "A total variation spectral framework for scale and texture analysis." SIAM Journal on Imaging Sciences 7.4 (2014): 1937-1961.

f = zeros(params.Method.n);
% CheckUps
if ~exist('params') || (exist('params') && ~isfield(params,'Method'))
    params.Method.Method_Type = 'AHMOD_TV';
    params.Method.Presets_Enable = 0; 
    params.Method.maxIter = 500;
    params.Method.print = 0;
elseif exist('params') && isfield(params,'Method')
    if  ~isfield(params.Method,'Method_Type'), params.Method.Method_Type = 'AHMOD_TV'; end
    if  ~isfield(params.Method,'Presets_Enable'), params.Method.Presets_Enable = 0; end
    if  ~isfield(params.Method,'maxIter'), params.Method.maxIter = 500; end
    if  ~isfield(params.Method,'print'), params.Method.print = 0; end
    if  ~isfield(params.Method,'k'), params.Method.k = 1; end
    if  ~isfield(params.Method,'CalcAInside'), params.Method.CalcAInside = 0; end
end
if ~exist('params') || (exist('params') && ~isfield(params,'ScaleParams'))
    params.ScaleParams.Max_time = 10;
    params.ScaleParams.Num_of_bands = 10;
    params.ScaleParams.dt = 1;
elseif exist('params') && isfield(params,'ScaleParams')    
    if  ~isfield(params.ScaleParams,'Max_time'), params.ScaleParams.Max_time = 10; end
    if  ~isfield(params.ScaleParams,'Num_of_bands'), params.ScaleParams.Num_of_bands = 10; end
    if  ~isfield(params.ScaleParams,'dt'), params.ScaleParams.dt = 1; end
end
if ~exist('params') || (exist('params') && ~isfield(params,'CalcAInside')), params.CalcAInside = 0; end


lambda = 1/(params.ScaleParams.dt);
NumIter = round(params.ScaleParams.Max_time/params.ScaleParams.dt);
Phi = zeros(params.Method.n);
Method_Type = params.Method.Method_Type;
A_Model = params.A_Model;
dt = params.ScaleParams.dt;
bi = b;
if params.Method.Presets_Enable
    if strcmp(Method_Type,'AHMOD_TV')
        [u0,lastErr(1),postsets{1}] = AHMOD_TV(A_Model,bi,lambda,params.Method); bi = A_Model*u0(:);
        [u1,lastErr(2),postsets{2}] = AHMOD_TV(A_Model,bi,lambda,params.Method,postsets{1}); bi = A_Model*u1(:);
        [u2,lastErr(3),postsets{3}] = AHMOD_TV(A_Model,bi,lambda,params.Method,postsets{2}); bi = A_Model*u2(:);
    elseif strcmp(Method_Type,'AHMOD_AATV')
        [u0,lastErr(1),postsets{1}] = AHMOD_TV(A_Model,bi,lambda,params.Method); bi = A_Model*u0(:);
        if ~isfield(params.Method,'A'), params.Method.A = CalcA( u0 ,params.Method.k); end
        [u1,lastErr(2),postsets{2}] = AHMOD_AATV(A_Model,bi,lambda,params.Method,postsets{1}); bi = A_Model*u1(:); 
        [u2,lastErr(3),postsets{3}] = AHMOD_AATV(A_Model,bi,lambda,params.Method,postsets{2}); bi = A_Model*u2(:);  
    end
else
    if strcmp(Method_Type,'AHMOD_TV')
        [u0,lastErr(1),postsets{1}] = AHMOD_TV(A_Model,bi,lambda,params.Method); bi = A_Model*u0(:);
        [u1,lastErr(2),postsets{2}] = AHMOD_TV(A_Model,bi,lambda,params.Method); bi = A_Model*u1(:);
        [u2,lastErr(2),postsets{3}] = AHMOD_TV(A_Model,bi,lambda,params.Method); bi = A_Model*u2(:);
    elseif strcmp(Method_Type,'AHMOD_AATV')
        [u0,lastErr(1),postsets{1}] = AHMOD_TV(A_Model,bi,lambda,params.Method); bi = A_Model*u0(:);
        if ~isfield(params.Method,'A'), params.Method.A = CalcA( u0 ,params.Method.k); end
        [u1,lastErr(2),postsets{2}] = AHMOD_AATV(A_Model,bi,lambda,params.Method); bi = A_Model*u1(:); 
        [u2,lastErr(3),postsets{3}] = AHMOD_AATV(A_Model,bi,lambda,params.Method); bi = A_Model*u2(:); 
    end
end
   
u(:,:,1) = u0;
u(:,:,2) = u1;
u(:,:,3) = u2;

% setting waitbar
h = waitbar(0,'Please wait...','Name','Approximating TV Flow...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

for i=1:NumIter
    startCyc = tic;
    ddu = (u0+u2-2*u1)/(dt*dt);  % one/two more iter
    t = i*dt;
    phi = ddu*t;
    Phi(:,:,i) = phi;
    S(i) = sum(abs(phi(:)));
    if (i<NumIter) % not last iteration
        u0=u1;
        u1=u2;
        if params.Method.Presets_Enable
            if strcmp(Method_Type,'AHMOD_TV')
                 bi = A_Model*u2(:); [u2,lastErr(i+3),postsets{i+3}] = AHMOD_TV(A_Model,bi,lambda,params.Method,postsets{i+1});
            elseif strcmp(Method_Type,'AHMOD_AATV')
                 bi = A_Model*u2(:); [u2,lastErr(i+3),postsets{i+3}] = AHMOD_AATV(A_Model,bi,lambda,params.Method,postsets{i+1});
            end
        else
            if strcmp(Method_Type,'AHMOD_TV')
                bi = A_Model*u2(:); [u2,lastErr(i+3),postsets{i+3}] = AHMOD_TV(A_Model,bi,lambda,params.Method); 
            elseif strcmp(Method_Type,'AHMOD_AATV')
                bi = A_Model*u2(:); [u2,lastErr(i+3),postsets{i+3}] = AHMOD_AATV(A_Model,bi,lambda,params.Method); 
            end
        end
        u(:,:,i+3) = u2;
        if params.CalcAInside && strcmp(Method_Type,'AHMOD_AATV')
            params.Method.A = CalcA( u2 ,params.Method.k);
        end
    end
    
    if getappdata(h,'canceling')
        NumIter = i;
        break
    end
    TCyc = toc(startCyc);
    rmnTotalTime = (NumIter - i)*TCyc;
    Hours = fix(rmnTotalTime/3600);
    Min = fix(mod(rmnTotalTime,3600)/60);
    Sec = fix(mod(rmnTotalTime,60));
    remainingTimeStr = [num2str(Hours) ':' num2str(Min) ':' num2str(Sec)];
    figure(100); plot((1:length(S))*dt,S);
    figure(101); imagesc(u2);colorbar;
    waitbar(i/NumIter,h,['Iter ' num2str(i) ' Out of ' num2str(NumIter) '. Remaning Time: ' remainingTimeStr])
end 
delete(h) 

f_r = (NumIter+1)*u1-NumIter*u2;  % residual image
T = (1:NumIter)*dt;

for j=1:NumIter
    p(:,:,j) = (u(:,:,j)-u(:,:,j+1))/dt;
end


XTV.S = S; XTV.T = T; XTV.Phi = Phi; XTV.p = p; XTV.u = u; XTV.f_r = f_r; XTV.b = b; XTV.lastErr = lastErr; 
XTV.Method = params.Method; XTV.ScaleParams = params.ScaleParams; XTV.postsets = postsets;
if  strcmp(Method_Type,'AHMOD_AATV'), XTV.A = params.Method.A; end