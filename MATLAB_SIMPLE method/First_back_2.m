function A = First_back_2(n,h)
v=ones(1,n-1);
v1=-1*ones(1,n-2);
A=diag(v,0)+diag(v1,-1);
v2=[zeros(1,n-2),-1];
A=[A;v2]/h;