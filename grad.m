function [fx,fy] = grad(P)
	fx = P(:,[2:end end])-P;
	fy = P([2:end end],:)-P;
end