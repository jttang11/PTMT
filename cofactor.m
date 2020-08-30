OrthorhombicTwinningEq;
%[m,l]=[1,2] for compound;
C1=trace(U1'*U1)-1/4*det(U1'*U1)-2-1/4*a1'*a1;
C2=a1'*(U1*(det(U1'*U1).*inv(U1'*U1-eye(3)))*n1);

OrthorhombicTwinningEq;
%[m,l]=[1,3]
Type1C1=trace(U1'*U1)-1/4*det(U1'*U1)-2-1/4*a1'*a1;
Type1C2=a1'*(U1*(det(U1'*U1).*inv(U1'*U1-eye(3)))*n1);
Type2C1=trace(U1'*U1)-1/4*det(U1'*U1)-2-1/4*a2'*a2;
Type2C2=a2'*(U1*(det(U1'*U1).*inv(U1'*U1-eye(3)))*n2);