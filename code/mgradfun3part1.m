function grad1=mgradfun3part1(pars,F,FF,data,Lambda0,d1,n,ranvec,nseed)
s=pars(1);
l=pars(2);
l0=pars(3);
L=l*Lambda0+l0*ones(n,1);

xqx=@(x) s*FF*x+myidct2(L.*mydct2(x,n,1),n,1);
b=s*(F')*data;
m=symmlq(xqx,b,1e-12,300);

%first derivtive of l with respect to s
l1_1=(1/2)*d1*s^(-1);%the first part of the equation
l1_3=(1/2)*(data-F*m)'*(data-F*m);%the third part of the equation

temp=zeros(nseed,1);
parfor i=1:1:nseed
    temp(i,1)=ranvec(:,i)'*symmlq(xqx,FF*ranvec(:,i),1e-10,200);%FF is dA/ds
end
l1_2=(1/2)*mean(temp);%the second part of the equation

l1=l1_1-l1_2-l1_3;

%first derivative of l with respect to l
xdx=@(x) L.*x;

temp=zeros(nseed,1);
parfor i=1:1:nseed
    temp(i,1)=ranvec(:,i)'*symmlq(xdx,Lambda0.*ranvec(:,i),1e-10,200);
end
l2_1=(1/2)*mean(temp);%the first part of the equation

temp=zeros(nseed,1);
parfor i=1:1:nseed
    temp(i,1)=ranvec(:,i)'*symmlq(xqx,myidct2(Lambda0.*mydct2(ranvec(:,i),n,1),n,1),1e-10,200);
end
l2_2=(1/2)*mean(temp);%the second part of the equation
l2_3=(1/2)*(mydct2(m,n,1)')*(Lambda0.*mydct2(m,n,1));%the third part of the equation
l2=l2_1-l2_2-l2_3;

%first derivative of l with respect to l0

temp=zeros(nseed,1);
parfor i=1:1:nseed
    temp(i,1)=ranvec(:,i)'*symmlq(xdx,ranvec(:,i),1e-10,200);
end
l3_1=(1/2)*mean(temp);

temp=zeros(nseed,1);
parfor i=1:1:nseed
   temp(i,1)=ranvec(:,i)'*symmlq(xqx,myidct2(mydct2(ranvec(:,i),n,1),n,1),1e-10,200);
end
l3_2=(1/2)*mean(temp);

l3_3=(1/2)*(mydct2(m,n,1)')*(mydct2(m,n,1));
l3=l3_1-l3_2-l3_3;

grad1=[l1;l2;l3];
end
