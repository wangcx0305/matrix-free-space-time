function iM=iPx(xx,n,t)
x=reshape(xx,n,t);
iMM=zeros(n,t);
for i=1:1:t
    iMM(:,i)=myidct2(x(:,i),n,1);
end
iM=reshape(iMM,n*t,1);
end