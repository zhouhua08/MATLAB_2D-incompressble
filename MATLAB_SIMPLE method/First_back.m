function A = First_back(n,h)
v=-1*ones(1,n-1);
v1=ones(1,n-2);
A=diag(v,0)+diag(v1,1);
v2=[zeros(n-2,1);1];
A=[A,v2]/h;