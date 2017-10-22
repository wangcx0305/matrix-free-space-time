c=128;%column size
r=128;%row size
n=c*r;%c*r
t=10;%time
F=speye(n*t);
FF=F'*F;
ssample=1;
lsample=2;
l0 = 1;
Dr=(2*(1-cos(pi*(0:(r-1))/r)))';
Dc=(2*(1-cos(pi*(0:(c-1))/c)))';
Lsample=0.5*lsample*(repmat(Dr,c,1)+kron(Dc,ones(r,1)))+0.01*ones(n,1);
I=diag(ones(n,1));
L1sample=(1./Lsample).*(diag(I)-exp(-l0 * Lsample));
iL1sample=1./L1sample;

rng(19920305);
z=myidct2(sqrt(1./Lsample).*normrnd(0,1,n,1),n,1);%theta0
thetasample=zeros(n*t,1);
thetasample(1:n,1)=myidct2(exp(-l0 * Lsample/2).*mydct2(z,n,1),n,1)+myidct2(sqrt(L1sample).*normrnd(0,1,n,1),n,1);
for i=2:1:t
thetasample((n*(i-1)+1):n*i,1)=myidct2(exp(-l0 * Lsample/2).*mydct2(thetasample((n*(i-2)+1):n*(i-1),1),n,1),n,1)+...
                                myidct2(sqrt(L1sample).*normrnd(0,1,n,1),n,1);
end
sigmasample=sqrt(1/ssample)*normrnd(0,1,n*t,1);
datasample=F*thetasample+sigmasample;

rng(19920305);
w=(rand(n,t)>0.2);
idx=find(w==1);
mdata64=nan(n,t);
mdata64(idx)=datasample(idx);
csvwrite('mdata128_new.csv',mdata64);
