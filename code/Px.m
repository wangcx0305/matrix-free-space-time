function M=Px(xx,n,t)
x=reshape(xx,n,t);
MM=zeros(n,t);
parfor i=1:1:t
    MM(:,i)=mydct2(x(:,i),n,1);
end
M=reshape(MM,n*t,1);
end