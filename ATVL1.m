function x = ATVL1(x0,params)

x = x0;

% line search parameters
maxlsiter = params.lineSearchItnlim ;
gradToll = params.gradToll ;
alpha = params.lineSearchAlpha; 
beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;
t = 1;

% copmute g0  = grad(Phi(x))
g0 = wGradient(x,params);

dx = -g0;


% iterations
while(1)

    % Calculate Matrix B
    params.B = CalcB( params.L1'*x , params.k);
    
    % backtracking line-search
	[MML1tx, MML1tdx, DL1tx, DL1tdx] = preobjective(x, dx, params);
	f0 = objective(MML1tx, MML1tdx, DL1tx, DL1tdx,x,dx, 0, params);
	t = t0;
        [f1, ERRobj, RMSerr,ATV,L1]  =  objective(MML1tx, MML1tdx, DL1tx, DL1tdx,x,dx, t, params);

	
	lsiter = 0;

	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 && (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		[f1, ERRobj, RMSerr,ATV,L1]  =  objective(MML1tx, MML1tdx, DL1tx, DL1tdx,x,dx, t, params);
	end

	if lsiter == maxlsiter
		disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
		return;
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2
		t0 = t0 * beta;
	end 
	
	if lsiter<1
		t0 = t0 / beta;
	end

	x = (x + t*dx);
    
    if params.noneg
    x=params.L1*(max(params.L1'*x,0));
    end

	%--------- uncomment for debug purposes ------------------------	
	disp(sprintf('%d   , obj: %f, RMS: %f, L-S: %f, %d', k,f1,RMSerr,lsiter));

	%---------------------------------------------------------------
	
    %conjugate gradient calculation
    
	[g1,gradObj,gradL1,gradATV]= wGradient(x,params);
     bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+1e-50);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
    
    if params.debug == 1
        figure(100);imagesc(params.L1'*x);colorbar;title('f_n');
        figure(101);imagesc(params.L1'*gradObj);colorbar;title('gObj');
        figure(102);imagesc(params.L1Weight*(params.L1'*gradL1));colorbar;title('gL1');
        figure(103);imagesc(params.ATVWeight*(params.L1'*gradATV));colorbar;title('gATV');
        figure(104);imagesc(params.L1'*g1);colorbar;title('gTotal');
        drawnow;
    end
	if (k > params.Itnlim) || (norm(dx(:)) < gradToll) 
		break;
	end

end


return;


function [MML1tx, MML1tdx, DL1tx, DL1tdx] = preobjective(x, dx, params)

% precalculates transforms to make line search cheap
shw=params.L1'*x;

shw1=params.L1'*dx;

MML1tx = params.MM*(shw(:));
MML1tdx = params.MM*(shw1(:));

if params.ATVWeight
    DL1tx = params.ATV(shw,params.B);
    DL1tdx = params.ATV(shw1,params.B);
   
else
    DL1tx = 0;
    DL1tdx = 0;
end





function [res, obj, RMS,ATV,L1] = objective(MML1tx, MML1tdx, DL1tx, DL1tdx, x,dx,t, params)
%calculated the objective function

p = params.pNorm;

obj = MML1tx + t*MML1tdx - params.data;
obj = obj(:)'*obj(:);

if params.ATVWeight
    w = DL1tx(:) + t*DL1tdx(:);
    ATV = (w.*conj(w)+params.l1Smooth).^(p/2); 
else
    ATV = 0;
end

if params.L1Weight
   w = x(:) + t*dx(:); 
   L1 = (w.*conj(w)+params.l1Smooth).^(p/2);
else
    L1=0;
end



ATV = sum(ATV.*params.ATVWeight(:));
L1 = sum(L1.*params.L1Weight(:));
RMS = sqrt(obj/sum(abs(params.data(:))>0));

res = obj + (ATV) + (L1) ;

function [grad,gradObj,gradL1,gradATV] = wGradient(x,params)

gradL1 = 0;
gradATV = 0;

gradObj = gOBJ(x,params);
if params.L1Weight
gradL1 = gL1(x,params);
end
if params.ATVWeight
gradATV = gATV(x,params);
end


grad = (gradObj +  params.L1Weight.*gradL1 + params.ATVWeight.*gradATV);




function gradObj = gOBJ(x,params)
% computes the gradient of the data consistency

	
    shw=params.L1'*x;
    shw1=params.MM'*(params.MM*(shw(:)) - params.data);
        gradObj = params.L1*(reshape(shw1,params.N,params.N));
       

%gradObj = 2*gradObj ;


function grad = gL1(x,params)
% compute gradient of the L1 transform operator

p = params.pNorm;

grad = p*x.*(x.*conj(x)+params.l1Smooth).^(p/2-1);


function grad = gATV(x,params)
% compute gradient of TV operator

p = params.pNorm;

Dx =params.ATV(params.L1'*x,params.B);
% G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);
G = p*Dx.*(sum(abs(Dx).^2,3) + params.l1Smooth).^(p/2-1);
grad = params.L1*(params.GATV(G,params.B));






