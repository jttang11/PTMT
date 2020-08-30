%generate the compostion versus lambda 2 matrix, named LAM
%expected input parameter lattpara=[row 1 concentration of alloying elements,
%row 2 lattice parameter of B2 phase, row 3-5 lattice constant of B19 phase/ row 3-6 lattice constant of B19' M-II martensite]
%eg [5,3.02,2.87,4.45,4.52] B2-orthohombic B19;
%[5,3.02,2.87,4.45,4.52,96.52] B2-Monoclinic II B19'
function [LAM]=lambdagenerator(lattpara)
m=size(lattpara);
for i = 1:m(1)
    LAM(i,1)=lattpara(i,1);
end
if (m(2))==5 %B2-B19 transformation
    alpha=lattpara(:,4)./(lattpara(:,2).*2^0.5);
    beta=lattpara(:,3)./lattpara(:,2);
    gamma=lattpara(:,5)./(lattpara(:,2).*2^0.5);
    for j = 1:m(1)
        [~,lambda2]=orthohombic(alpha(j),beta(j),gamma(j));
        LAM(j,2)=lambda2;
    end
end
if (m(2))==6 %B2-B19' Monoclinic II transformation
    p1=lattpara(:,3)./lattpara(:,2);
    p2=lattpara(:,4)./lattpara(:,2);
    p3=lattpara(:,5)./lattpara(:,2);
    p4=-cosd(lattpara(:,6));  %cos beta angle in degrees
    for j = 1:m(1)
        [~,lambda2]=monoclinic2(p1(j),p2(j),p3(j),p4(j));
        LAM(j,2)=lambda2;
    end
end
    