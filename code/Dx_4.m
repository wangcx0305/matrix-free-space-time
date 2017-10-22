function C = Dx_4(uu, t, n, p0, L, iL1)
%uu is the mean vector of theta
%t is for time indicator, n is for the number of observation
%L and iL1 is Lambda and Lambda1^{-1}
u = reshape(uu, n, t);
D1 = @(x) iL1 .* x;%Lambda1^{-1}x
D2 = @(x) -exp(-(p0 * L) / 2) .* iL1 .* x;%-e^{-Lambda/2}*Lambda1^{-1}*x
D3 = @(x) exp(-(p0 * L)) .* iL1 .* x;%e^{-Lambda}*Lambda1^{-1}*x
%C(1:n,1)=D1(u(:,1))+D2(u(:,2));
C(1:n, 1) = D1(u(:, 1)) + D2(u(:, 2));
for i = 2:1:(t-1)
    C((n*(i-1)+1):n*i,1) = D2(u(:,i-1)) + D1(u(:,i)) + D3(u(:,i)) + D2(u(:,i+1));
end
C((n*(t-1)+1):n*t,1) = D2(u(:,t-1)) + D1(u(:,t));
end