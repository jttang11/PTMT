%Solve the compatiability equation
%Must have run twinning equation mfile first
%Firstly check whether there exists any solution
%First pair of solution
delta1=a1'*(G*inv(G*G-eye(3))*n1);
eta1=trace(G*G)-det(G*G)-2+(a1'*a1)/(2*delta1);
if (delta1<-2 || eta1>0)
    lambda=1/2*(1-sqrt(1+2/delta1));
    C=(G+lambda*n1*a1')*(G+lambda*a1*n1');
    %ETG records the eigen values of C matrix
    [vec,EIG]=eig(C);
    lambda1=EIG(1,1);
    lambda2=EIG(2,2);
    lambda3=EIG(3,3);
    e1=vec(:,1);
    e2=vec(:,2);
    e3=vec(:,3);
    %Then the solution of compatiability equation
    M1=(sqrt(lambda3)-sqrt(lambda1))/sqrt(lambda3-lambda1).*(-sqrt(1-lambda1).*e1+sqrt(lambda3-1).*e3);
    rho1=norm(M1);
    m1=M1/rho1;
    b1=rho1*(sqrt(lambda3*(1-lambda1)/(lambda3-lambda1))*e1+sqrt(lambda1*(lambda3-1)/(lambda3-lambda1))*e3);
    %2nd solution
    M2=(sqrt(lambda3)-sqrt(lambda1))/sqrt(lambda3-lambda1).*(-sqrt(1-lambda1).*e1-sqrt(lambda3-1).*e3);
    rho2=norm(M2);
    m2=M2/rho2;
    b2=rho2*(sqrt(lambda3*(1-lambda1)/(lambda3-lambda1))*e1-sqrt(lambda1*(lambda3-1)/(lambda3-lambda1))*e3);
else
    fprintf('The variant pair is not compatible with the parent phase\n');
end
%%
if (delta1<-2 || eta1>0)
    fprintf('For a1 and n1\n');
    fprintf('The volume fraction of 1 variant is %.4f\n',lambda);
    fprintf('m1=[%.3f,%.3f,%.3f]\t',m1(1),m1(2),m1(3));
    fprintf('b1=[%.3f,%.3f,%.3f]\n',b1(1),b1(2),b1(3));
    fprintf('m2=[%.3f,%.3f,%.3f]\t',m2(1),m2(2),m2(3));
    fprintf('b2=[%.3f,%.3f,%.3f]\n',b2(1),b2(2),b2(3));
end
%%
delta2=a2'*(G*inv(G*G-eye(3))*n2);
eta2=trace(G*G)-det(G*G)-2+(a2'*a2)/(2*delta2);
if (delta2<-2 || eta2>0)
    lambda=1/2*(1-sqrt(1+2/delta2));
    C=(G+lambda*n2*a2')*(G+lambda*a2*n2');
    %ETG records the eigen values of C matrix
    [vec,EIG]=eig(C);
    lambda1=EIG(1,1);
    lambda2=EIG(2,2);
    lambda3=EIG(3,3);
    e1=vec(:,1);
    e2=vec(:,2);
    e3=vec(:,3);
    %Then the solution of compatiability equation
    M3=(sqrt(lambda3)-sqrt(lambda1))/sqrt(lambda3-lambda1).*(-sqrt(1-lambda1).*e1+sqrt(lambda3-1).*e3);
    rho3=norm(M1);
    m3=M3/rho3;
    b3=rho3*(sqrt(lambda3*(1-lambda1)/(lambda3-lambda1))*e1+sqrt(lambda1*(lambda3-1)/(lambda3-lambda1))*e3);
    %2nd solution
    M4=(sqrt(lambda3)-sqrt(lambda1))/sqrt(lambda3-lambda1).*(-sqrt(1-lambda1).*e1-sqrt(lambda3-1).*e3);
    rho4=norm(M4);
    m4=M4/rho4;
    b4=rho4*(sqrt(lambda3*(1-lambda1)/(lambda3-lambda1))*e1-sqrt(lambda1*(lambda3-1)/(lambda3-lambda1))*e3);
else
    fprintf('The variant pair is not compatible with the parent phase\n');
end
%%
if (delta2<-2 || eta2>0)
    fprintf('For a2 and n2\n');
    fprintf('The volume fraction of 1 variant is %.4f\n',lambda);
    fprintf('m3=[%.3f,%.3f,%.3f]\t',m3(1),m3(2),m3(3));
    fprintf('b3=[%.3f,%.3f,%.3f]\n',b3(1),b3(2),b3(3));
    fprintf('m4=[%.3f,%.3f,%.3f]\t',m4(1),m4(2),m4(3));
    fprintf('b4=[%.3f,%.3f,%.3f]\n',b4(1),b4(2),b4(3));
end
