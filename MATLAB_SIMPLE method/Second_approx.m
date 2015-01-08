% K1_full stores the full matrix with a tridiagonal structure
% approximates the second order derivative
function A = K1_full(n,h)
v=-2*ones(1,n);
A=diag(v,0)/h^2;
