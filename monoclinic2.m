function [U, lambda2]=monoclinic2(p1,p2,p3,p4)
%monoclinic would calculate the B19' phase eigen value
%Monoclic-II martensite corresponds to <100>cubic
%p1=a/a0,p2=b/a0,p3=c/a0;p4=-cos(β)
%U is the transformation stretch tensor
%lambda2 is the 2nd eigenvalue
U=[p1,0,(p1*p3*p4)^0.5;
   0,p2/sqrt(2),0;
   (p1*p3*p4)^0.5,0,p3*sqrt(1+p4^2/2)];
%x denote the matrix of eigen vectors while y denotes the eigen matrix
[~,y]=eig(U);
lambda2=y(2,2);
