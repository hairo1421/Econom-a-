%DGP AR(1)
clear;
T=200;
c=1;
phi1=0.5;
s2=1;
y(1,1)=random('normal',c/(1-phi1),sqrt(s2/(1-phi1^2)));
for t=2:T
    y(t,1)=c+phi1*y(t-1,1)+randn;
end
%OLS
X=[ones(T-1,1),y(1:end-1,:)];
Y=y(2:end,:);
bols=inv(X'*X)*X'*Y;
%MCMC
%prior
beta0=zeros(2,1);
sigma0=[10,0;0,0.5];
A0=1/2;
B0=2/2;
S2=2;
N0=2000;
for j=1:10000
    sigmapost=inv(inv(sigma0)+inv(S2*inv(X'*X)));
    meanpost=sigmapost*(inv(sigma0)*beta0+inv(S2*inv(X'*X))*bols);
    C=chol(sigmapost)';
    beta=meanpost+C*randn(2,1);
    SEC=sum((Y-X*beta).^2);
    A=A0+SEC;
    B=B0+T-1-1;
    S2=gamrnd(A,B^-1)^-1;
    if j>N0
        betadraw(:,j-N0)=beta;
        s2draw(:,j-N0)=S2;
    end
end

