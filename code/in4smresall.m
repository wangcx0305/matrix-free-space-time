smres = csvread('smres.csv');
smres = smres(:, 1:12);
c1 = 20;%column size
r1 = 40;%row size
c2 = 80;
r2 = 160;
n = c2 * r2;
t = 12;
temp = smres(:, 1);
F1 = zeros(length(temp), n);
for i = 1:1:length(temp)
    a = ceil(i / c1) - 1;
    b = a * 4 * c2 + (i - a * c1 - 1) * 4 + 1;
    F1(i, b) = 1/16;
    F1(i, b + 1) = 1/16;
    F1(i, b + 2) = 1/16;
    F1(i, b + 3) = 1/16;
    F1(i, b + c2) = 1/16;
    F1(i, b + c2 + 1) = 1/16;
    F1(i, b + c2 + 2) = 1/16;
    F1(i, b + c2 + 3) = 1/16;
    F1(i, b + 2 * c2) = 1/16;
    F1(i, b + 2 * c2 + 1) = 1/16;
    F1(i, b + 2 * c2 + 2) = 1/16;
    F1(i, b + 2 * c2 + 3) = 1/16;
    F1(i, b + 3 * c2) = 1/16;
    F1(i, b + 3 * c2 + 1) = 1/16;
    F1(i, b + 3 * c2 + 2) = 1/16;
    F1(i, b + 3 * c2 + 3) = 1/16;
    
end
ord = find(isnan(temp));
F1(ord, :) = [];
F1 = sparse(F1);
F = sparse(kron(eye(t), F1));
FF = F' * F;
temp = reshape(smres, 800 * t, 1);
datasample = temp(~isnan(temp));
d1 = length(datasample);
Dr = (2 * (1 - cos(pi * (0:(r2 - 1)) / r2)))';
Dc = (2 * (1 - cos(pi * (0:(c2 - 1)) / c2)))';
Lambda0 = 0.5 * repmat(Dr,c2,1) + 0.5 * kron(Dc,ones(r2,1));

rng(19920305)
nseed = 50;
ranvec=2 * (rand(n * t, nseed) < 0.5) - 1;

options=optimset('Display','iter','Tolfun',0.01,'TolX',0.001,'MaxFunEvals',500,'MaxIter',300);
pars0=[1.5, 200, 50, 1];

gradiso=@(pars) mgradfun4par(pars,F,FF,datasample,Lambda0,d1,n,t,ranvec,nseed);
[x,fval,exitflag,output,hess]=fsolve(gradiso,pars0,options);
x(1)
x(2)
x(3)
x(4)
in = -hess
v = inv(-hess);
se = sqrt(diag(v))
se(1) 
se(2)
se(3)
se(4)

[V, D] = eig(in);
Dnew = zeros(4,4);
for i = 1:1:4
     if real(D(i, i)) > 0
         Dnew(i, i) = real(D(i, i));
     end
         
end
in = V * Dnew * V';
v = pinv(in)
se = sqrt(diag(v))
se(1)
se(2)
se(3)
se(4)











