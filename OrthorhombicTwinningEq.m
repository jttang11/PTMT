%compute the orthorhombic twinning equation solutions
%a vector is the direction of shear 
%n vector could be related to the shear plane by Uj_inverse
%e vector is the rotation axis
%James version

%Firstly input lattice parameters
a0=input('Please input the cubic phase lattice constant: ');
a=input('Please input the orthorhombic phase lattice constant a: ');
b=input('Please input the orthorhombic phase lattice constant b: ');
c=input('Please input the orthorhombic phase lattice constant c: ');

%Then define alpha, beta, gamma
alpha=a/(a0*2^0.5);
beta=b/a0;
gamma=c/(a0*2^0.5);

%Define Variant 1-6 (Bhattacharya P27 table4.3)
U1=[(alpha+gamma)/2,0,(alpha-gamma)/2;
    0,beta,0;
    (alpha-gamma)/2,0,(alpha+gamma)/2];
U2=[(alpha+gamma)/2,0,(gamma-alpha)/2;
    0,beta,0;
    (gamma-alpha)/2,0,(alpha+gamma)/2];
U3=[(alpha+gamma)/2,(alpha-gamma)/2,0;
    (alpha-gamma)/2,(alpha+gamma)/2,0;
    0,0,beta];
U4=[(alpha+gamma)/2,(gamma-alpha)/2,0;
    (gamma-alpha)/2,(alpha+gamma)/2,0;
    0,0,beta];
U5=[beta,0,0;
    0,(alpha+gamma)/2,(alpha-gamma)/2;
    0,(alpha-gamma)/2,(alpha+gamma)/2];
U6=[beta,0,0;
    0,(alpha+gamma)/2,(gamma-alpha)/2;
    0,(gamma-alpha)/2,(alpha+gamma)/2];
%%
%Input lattice variant pair, expected value [index of variant 1,index of variant 2]
variantpair=input('Please input vairant pair in [m,l]: ');
flag=0;
while flag==0
    if (ismember(variantpair(1),[1,2,3,4,5,6])...
        &&...
        ismember(variantpair(2),[1,2,3,4,5,6])...
        &&...
        variantpair(1)~=variantpair(2))
        flag=1;
        G=eval(strcat('U',num2str(variantpair(1))));
        F=eval(strcat('U',num2str(variantpair(2))));
    else
        disp('Input error')
        variantpair=input('Please input vairant pair in [m,l]: ');
    end
end
%%
%So far we have got G matrix and twinning pairs, we will now use e vectors
%to solve the twinning equation Bhattacharya note p35
%The output will be 2 set of solutions (n1,a1) and (n2,a2)   
if (all(variantpair==[1,2]) || all(variantpair==[2,1])) %variant pair 1,2
    e=transpose([1 0 0]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end

if (all(variantpair==[1,3]) || all(variantpair==[3,1])) %variant pair 1,3
    e=transpose(1/sqrt(2).*[0 1 1]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end

if (all(variantpair==[1,4]) || all(variantpair==[4,1])) %variant pair 1,4
    e=transpose(1/sqrt(2).*[0 -1 1]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end

if (all(variantpair==[1,5]) || all(variantpair==[5,1])) %variant pair 1,5
    e=transpose(1/sqrt(2).*[1 1 0]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end

if (all(variantpair==[1,6]) || all(variantpair==[6,1])) %variant pair 1,6
    e=transpose(1/sqrt(2).*[1 -1 0]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end

if (all(variantpair==[2,3]) || all(variantpair==[3,2])) %variant pair 2,3
    e=transpose(1/sqrt(2).*[0 -1 1]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end

if (all(variantpair==[2,4]) || all(variantpair==[4,2])) %variant pair 2,4
    e=transpose(1/sqrt(2).*[0 1 1]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end

if (all(variantpair==[2,5]) || all(variantpair==[5,2])) %variant pair 2,5
    e=transpose(1/sqrt(2).*[1 -1 0]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end

if (all(variantpair==[2,6]) || all(variantpair==[6,2])) %variant pair 2,6
    e=transpose(1/sqrt(2).*[1 1 0]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end

if (all(variantpair==[3,4]) || all(variantpair==[4,3])) %variant pair 3,4
    e=transpose([0 1 0]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end

if (all(variantpair==[3,5]) || all(variantpair==[5,3])) %variant pair 3,5
    e=transpose(1/sqrt(2).*[1 0 -1]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end

if (all(variantpair==[3,6]) || all(variantpair==[6,3])) %variant pair 3,6
    e=transpose(1/sqrt(2).*[1 0 1]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end

if (all(variantpair==[4,5]) || all(variantpair==[5,4])) %variant pair 4,5
    e=transpose(1/sqrt(2).*[1 0 -1]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end

if (all(variantpair==[4,6]) || all(variantpair==[6,4])) %variant pair 4,6
    e=transpose(1/sqrt(2).*[1 0 1]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end

if (all(variantpair==[5,6]) || all(variantpair==[6,5])) %variant pair 5,6
    e=transpose([0 0 1]);
    n1=e;
    GT=transpose(inv(G));
    a1=2.*(GT*e/(norm(GT*e))^2-G*e);
    %now work with solution 2
    sol2=2.*(e-G'*G*e/(norm(G*e))^2);
    rho=norm(sol2);
    n2=sol2./rho;
    a2=rho.*(G*e);
end
%%
fprintf('solution 1: a1=[%.3f %.3f %.3f]\t',a1(1),a1(2),a1(3))
fprintf('n1=[%.3f %.3f %.3f]\n',n1(1),n1(2),n1(3))
fprintf('solution 2: a2=[%.3f %.3f %.3f]\t',a2(1),a2(2),a2(3))
fprintf('n2=[%.3f %.3f %.3f]\n',n2(1),n2(2),n2(3))