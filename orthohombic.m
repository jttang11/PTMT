function [U, lambda2]=orthohombic(alpha,beta,gamma)
%orthohombic would calculate the B19 phase eigen value
%U is the transformation stretch tensor
%lambda2 is the 2nd eigenvalue
U=[beta,0,0;
   0,(alpha+gamma)/2,(alpha-gamma)/2;
   0,(alpha-gamma)/2,(alpha+gamma)/2];
%x denote the matrix of eigen vectors while y denotes the eigen matrix
[x,y]=eig(U);
lambda2=y(2,2);
