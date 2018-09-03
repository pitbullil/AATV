function res = init()
% param = init()

res.MM = []; % The measurement operator 
res.L1 = []; % Sparse transform operator
res.TV = []; 	% the Total variation operator
res.LAP = []; 	% the Laplacian operator
res.data = []; % measurements to reconstruct from

res.TVWeight = 0.01;	% TV penalty
res.LAPWeight = 0.01;	% LAP penalty
res.L1Weight = 0.01;   % transform l1 penalty

res.Itnlim = 150;	% default number of iterations
res.gradToll = 1e-10;	% step size tollerance stopping criterea 

res.l1Smooth = 1e-15;	% smoothing parameter of L1 norm
res.pNorm = 1;  % type of norm to use (i.e. L1 L2 etc)

% line search parameters
res.lineSearchItnlim = 450;
res.lineSearchAlpha = 0.4;
res.lineSearchBeta = 0.1;
% res.lineSearchAlpha = 0.01;
% res.lineSearchBeta = 0.6;
res.lineSearchT0 = 1 ; % step size to start with




