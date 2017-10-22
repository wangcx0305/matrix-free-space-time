function grad = mgradfun4par(pars, F, FF, data, Lambda0, d1, n, t, ranvec, nseed)
p0 = pars(1);
p1 = pars(2);
p2 = pars(3);
p3 = pars(4);
I = speye(n);
L = p1 * Lambda0 + p2 * ones(n,1);%Lambda
L1 = (1./L) .* (diag(I) - exp(-p0 * L));%Lambda1
iL1 = 1./L1;%Lambda1^{-1}
b = p3 * (F') * data;
B1 = @(v1) Ax_4(v1, t, n, p3, p0, FF, L, iL1);%calculate X'Sigma^{-1}XY when Y is unknown
m = symmlq(B1, b, 1e-12, 300);
     
%first derivtive of l with respect to s
l3_1 = (1/2) * d1 * p3^(-1);%the first part of the equation
l3_3 = (1/2) * (data - F*m)' * (data - F*m);%the third part of the equation
     
temp = zeros(nseed,1);
parfor i = 1:1:nseed
temp(i,1) = ranvec(:,i)' * symmlq(B1, FF * ranvec(:,i),1e-8,200);%FF is dA/ds
end
l3_2 = (1/2) * mean(temp);%the second part of the equation
l3 = l3_1 - l3_2 - l3_3;

%first derivative of l with respect to l
%first derivative of Lambda with respect to lambda2 which is l
D10 = -(iL1 .^ 2) .* exp(-p0 * L);%first derivative of Lambda1^{-1} with respect to lambda which is l
D20 = (1 + exp(-p0 * L)) .* D10 - exp(-p0 * L) .* iL1 .* L;%first derivative of Lambda1^{-1}+e^{-Lambda}*Lambda1^{-1} with respect to lambda which is l
D30 = exp(-p0 * L / 2) .* (D10 - (1/2) * iL1 .* L);%first derivative of e^{-Lambda/2}*Lambda_1^{-1} with respect to lambda which is l

G2=@(v2) Dx_4(v2, t, n, p0, L, iL1);

temp = zeros(nseed,1);
parfor i = 1:1:nseed
    temp(i,1) = ranvec(:,i)' * symmlq(G2, pDplx(ranvec(:,i), t, n, D10, D20, D30), 1e-8, 200);
end
l0_1 = (1/2) * mean(temp);%the first part of the equation

temp = zeros(nseed,1);
parfor i = 1:1:nseed
    c = iPx(pDplx(Px(ranvec(:,i),n,t),t,n,D10,D20,D30),n,t);
    temp(i,1)=(ranvec(:,i)')*symmlq(B1,c,1e-8,200);
end
l0_2 = (1/2) * mean(temp);%the second part of the equation

l0_3 = (1/2)*(Px(m,n,t)')*pDplx(Px(m,n,t),t,n,D10,D20,D30);%the third part of the equation

l0 = l0_1-l0_2-l0_3;


%first derivative of l with respect to l
D = Lambda0;%first derivative of Lambda with respect to lambda2 which is l
D11 = -(1 - exp(-p0 * L)).^(-2).*exp(-p0 * L).*(p0 * L).*D + (1./(1-exp(-p0 * L))).*D;%first derivative of Lambda1^{-1} with respect to lambda which is l
D21 = (1 + exp(-p0 * L)).*D11 - exp(-p0 * L).*iL1.* (p0 * D);%first derivative of Lambda1^{-1}+e^{-Lambda}*Lambda1^{-1} with respect to lambda which is l
D31 = exp(- p0 * L / 2).*(D11 - (p0 / 2).*iL1.*D);%first derivative of e^{-Lambda/2}*Lambda_1^{-1} with respect to lambda which is l

G2 = @(v2) Dx_4(v2,t,n,p0,L,iL1);

temp = zeros(nseed,1);
parfor i = 1:1:nseed
    temp(i,1) = ranvec(:,i)' * symmlq(G2,pDplx(ranvec(:,i),t,n,D11,D21,D31),1e-8,200);
end
l1_1 = (1/2) * mean(temp);%the first part of the equation

temp = zeros(nseed,1);
parfor i = 1:1:nseed
    c = iPx(pDplx(Px(ranvec(:,i),n,t),t,n,D11,D21,D31),n,t);
    temp(i,1) = (ranvec(:,i)') * symmlq(B1,c,1e-8,200);
end
l1_2 = (1/2) * mean(temp);%the second part of the equation

l1_3 = (1/2) * (Px(m,n,t)') * pDplx(Px(m,n,t),t,n,D11,D21,D31);%the third part of the equation

l1 = l1_1 - l1_2 - l1_3;

%first derivative of l with respect to l0
D12 = -(1 - exp(-p0 * L)).^(-2) .* exp(-p0 * L).*(p0 * L) + (1 ./ (1 - exp(-p0 * L)));%first derivative of Lambda1^{-1} with respect to lambda which is l
D22 = (1 + exp(-p0 * L)) .* D12 - p0 * exp(-p0 * L) .* iL1;%first derivative of Lambda1^{-1}+e^{-Lambda}*Lambda1^{-1} with respect to lambda which is l
D32 = exp(-p0 * L / 2) .* (D12 - (p0/2) * iL1);%first derivative of e^{-Lambda/2}*Lambda_1^{-1} with respect to lambda which is l

temp = zeros(nseed,1);
parfor i = 1:1:nseed
    temp(i,1) = ranvec(:,i)' * symmlq(G2,pDplx(ranvec(:,i),t,n,D12,D22,D32),1e-8,200);
end
l2_1 = (1/2) * mean(temp);%the first part of the equation

temp = zeros(nseed,1);
parfor i = 1:1:nseed
    c = iPx(pDplx(Px(ranvec(:,i),n,t),t,n,D12,D22,D32),n,t);
    temp(i,1) = (ranvec(:,i)')*symmlq(B1,c,1e-8,200);
end
l2_2 = (1/2) * mean(temp);%the second part of the equation

l2_3 = (1/2) * (Px(m,n,t)') * pDplx(Px(m,n,t),t,n,D12,D22,D32);%the third part of the equation

l2 = l2_1 - l2_2 - l2_3;

grad = [l3, l0, l1, l2];
end