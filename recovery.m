m=size(LoadD);
len=m(1);
i=1;
result=zeros(len,2);
while i < len
    e=LoadD(i,:);
    e=e./norm(e);
    t1=sqrt(e*(U1*e'));
    t2=sqrt(e*(U2*e'));
    t3=sqrt(e*(U3*e'));
    t4=sqrt(e*(U4*e'));
    t5=sqrt(e*(U5*e'));
    t6=sqrt(e*(U6*e'));
    result(i,1)=max([t1,t2,t3,t4,t5,t6])-1;
    result(i,2)=abs(min([t1,t2,t3,t4,t5,t6])-1);
    i=i+1;
end