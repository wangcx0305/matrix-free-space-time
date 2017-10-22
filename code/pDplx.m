%function to compute dD/d(lambda)*x, where D is the second part of
%covariance matrix, lambda could be lambda_2,lambda_3 or lambda_4
function PC=pDplx (uu,t,n,D1,D2,D3)%
%uu is the mean vector of theta
%t is for time indicator, n is for the number of observation
%D1,D2 and D3 are from first derivative of D with respect to lambda2
u = reshape(uu,n,t);
pDpl1=@(x) D1.*x;%D1*x
pDpl2=@(x) D2.*x;%D2*x
pDpl3=@(x) -D3.*x;%?D3*x
PC(1:n,1)=pDpl1(u(:,1))+pDpl3(u(:,2));
for i=2:1:(t-1)
    PC((n*(i-1)+1):n*i,1)=pDpl3(u(:,i-1))+pDpl2(u(:,i))+pDpl3(u(:,i+1));
end
PC((n*(t-1)+1):n*t,1)=pDpl3(u(:,t-1))+pDpl1(u(:,t));
end