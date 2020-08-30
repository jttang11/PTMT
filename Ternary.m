%Generate the lambda2 of a series of alloys
X=(40:0.1:60); %atomic ratio of Ni
Z=(0:0.1:20); %atomic ration of Pd
i=1;
sizeX=size(X);
sizeZ=size(Z);
Grid=zeros(40401,3);
%mesh the grid
while i <=sizeX(2)
    j=1;
    while j <= sizeZ(2)
        Grid(sizeZ(2)*(i-1)+j,1)=X(i);
        Grid(sizeZ(2)*(i-1)+j,3)=Z(j);
        Grid(sizeZ(2)*(i-1)+j,2)=100-X(i)-Z(j);
        j=j+1;
    end
    i=i+1;
end
gridsize=size(Grid);
lengrid=gridsize(1);
%%
%generate lattice parameter vector
k=1;
LatticeParameter=zeros(40401,5);
while k<lengrid
    LatticeParameter(k,1)=k;
    LatticeParameter(k,2)=3.018-0.3524*(Grid(k,1)-50)-0.3552*(Grid(k,2)-50)-0.3624*Grid(k,3); %a0
    LatticeParameter(k,3)=2.889-0.00282*Grid(k,3); %a
    LatticeParameter(k,4)=4.12-3.3192*(Grid(k,1)-50)-3.3592*(Grid(k,2)-50)-3.2642*Grid(k,3); %b
    LatticeParameter(k,5)=4.5269-0.4320*(Grid(k,1)-50)-0.4194*(Grid(k,2)-50)-0.4256*Grid(k,3); %c
    k=k+1;
end
%%
%now generate lambda2
[LAM]=lambdagenerator(LatticeParameter);