mdatat10 = csvread('mdata128_new.csv');
c = 128;%column size
r = 128;%row size
n = 128 * 128;%c*r
temp = mdatat10(:, 1);
idx = find(temp>0|temp<0);
d1 = length(idx);
F = sparse(1:d1,idx,1,d1,n);
datasample = temp(temp>0|temp<0);
FF = F'*F;
Dr = (2*(1-cos(pi*(0:(r-1))/r)))';
Dc = (2*(1-cos(pi*(0:(c-1))/c)))';
Lambda0 = 0.5*repmat(Dr,c,1) + 0.5*kron(Dc,ones(r,1));

rng(19920305)
nseed = 50;
ranvec = 2*(rand(n,nseed)<0.5)-1;

options = optimset('Display','iter','Tolfun',0.01,'TolX',0.001,'MaxFunEvals',500,'MaxIter',100);
pars0 = [1,2,0.01];

gradiso = @(pars) mgradfun3part1(pars,F,FF,datasample,Lambda0,d1,n,ranvec,nseed);
[x,fval,exitflag,output,hess] = fsolve(gradiso,pars0,options);
1 / x(1)
x(1)
x(2)
x(3)
v = inv(-hess)
se = sqrt(diag(v))
se(1)
se(1)/((x(1))^2)
se(2)
se(3)