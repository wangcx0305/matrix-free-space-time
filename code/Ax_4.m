function B = Ax_4(uu, t, n, p3, p0, FF, L, iL1)
%uu is the mean vector of theta
%t is for time indicator, n is for the number of observation
%l2 is sigma^{2}
%F is the observation matrix with dimension which is composed by
%F_1,...F_m, each F_i is n*n
%L and iL1 is Lambda and Lambda1^{-1}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = reshape(uu, n, t);
A1 = @(x) myidct2(iL1 .* mydct2(x, n, 1), n, 1);%V^{-1}x
A2 = @(x) -myidct2(exp(-(p0 * L) / 2) .* iL1 .* mydct2(x, n, 1), n, 1);%-G^TV^{-1}x
A3 = @(x) myidct2(exp(-(p0 * L) / 2) .* iL1 .* exp(-(p0 * L) / 2) .* mydct2(x,n,1),n,1);%G^{T}V^{-1}Gx
B1 = p3 * FF * uu;%sigma^{-2}F'F\beta
%B2(1:n,1)=A1(u(:,1))+A2(u(:,2));
B2(1:n, 1) = A1(u(:, 1)) + A2(u(:, 2));
for i = 2:1:(t-1)
    B2((n*(i-1)+1):n*i, 1) = A2(u(:,i-1)) + A1(u(:,i)) + A3(u(:,i)) + A2(u(:,i+1));
end
B2((n*(t-1)+1):n*t, 1) = A2(u(:,t-1)) + A1(u(:,t));
B = B1 + B2;
end