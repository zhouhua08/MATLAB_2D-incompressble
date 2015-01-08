% K1_full stores the full matrix with a tridiagonal structure
% represents the second order derivative
function A = K1_full(n,h,a11,a22)
% a11: Neumann=1, Dirichlet=2, Dirichlet mid=3;
v=[-1*a11,-2*ones(1,n-2),-1*a22];
v1=1*ones(1,n-1);
A=(diag(v,0)+diag(v1,1)+diag(v1,-1))/h^2;
