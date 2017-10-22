smres = csvread('smres.csv');
smres = smres(:, 1:12);
d = 1/12;
[n, t] = size(smres);
idx = find(~isnan(smres));
d1 = length(idx);
F = sparse(1:d1,idx,1,d1,n*t);
datasample = smres(~isnan(smres));
c = 20;%column size
r = 40;%row size
FF = F' * F;
Dr=(2 * (1 - cos(pi * (0:(r-1)) / r)))';
Dc=(2 * (1 - cos(pi * (0:(c-1)) / c)))';
Lambda0 = 0.5 * (repmat(Dr,c,1) + kron(Dc,ones(r,1)));

rng(19920305)
nseed = 50;
ranvec=2 * (rand(n * t, nseed) < 0.5) - 1;

options=optimset('Display','iter','Tolfun',0.01,'TolX',0.001,'MaxFunEvals',1000,'MaxIter',300);
pars0=[1.5, 100, 100, 1];

gradiso=@(pars) mgradfun4par(pars,F,FF,datasample,Lambda0,d1,n,t,ranvec,nseed);
[x,fval,exitflag,output,hess]=fsolve(gradiso,pars0,options);
x(1)
x(2)
x(3)
x(4)
in = -hess;
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






