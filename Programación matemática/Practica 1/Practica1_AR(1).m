clear
%DGP AR(1)
T=100;
c=1.5;
b1=.6;
sigerr=2;
%Momentos
mu=c/(1-b1);
sigma=sigerr/(1-b1^2);
%creando observaciones
%Nota:Debes ponerlo transpuesto
y=normrnd(mu,sigma);
%Generando la serie
for t=2:T
    y(t)=c+b1*y(t-1)+rand;
end 
 plot(y);
%comprpobra
X=[ones(T-1,1),y(1:end-1)'];
Y=y(2:end)';
ols=inv(X'*X)*X'*Y
 %MCMC
 
 %Momentos Prior
 mu0=[0;0];
 sig0=[15,0;0,.5];
 A0=1/2
 B0=2/2
 %Momentos likelihood
mulik=ols;
siglik=3.2;

for j=1:8000
 sigpos=inv(inv(sig0)+inv(siglik));
 mupos=sigpos*(inv(sig0)*mu0+inv(siglik)*mulik);
 theta(j,:)=mvnrnd(mupos,sigpos);
 ssr=sum((Y-X*ols).^2);
 A=A0+ssr;
 B=B0+T-2;
 s2draw(j,1)=gamrnd(A,B^-1)^-1;
 s2error=s2draw(j,1);
end 
n0=500;
%plot
subplot(2,2,1);
ksdensity(s2draw(n0+1:end,1));
title('varianza del error posterior');
subplot(2,2,2);
ksdensity(theta(n0+1:end,1));
title('intercepto posterior');
subplot(2,2,3);
ksdensity(theta(n0+1:end,2));
title('autoregresivo 1 posterior');
%HPDI
intsup=mupos(2)+1.96*sqrt(sigpos(4));
intinf=mupos(2)-1.96*sqrt(sigpos(4));
int=[intinf,mupos(2),intsup];
intsup_c=mupos(1)+1.96*sqrt(sigpos(1));
intinf_c=mupos(1)-1.96*sqrt(sigpos(1));
int_c=[intinf_c,mupos(1),intsup_c];


prior=mvnpdf(theta,mu0',sig0);
plot(prior)
title('prior');
prior2=gampdf(s2draw,A0,B0);
siglik2=[siglik,0;0,3];
plot(prior2)
title('Varianza prior')

lik=0
for j=1:8000
lik(j,:)=sum(log(mvnpdf(theta(j,:),mulik',siglik2)));
end
plot(lik)
title('likelihood')
[cc,rr] = max(lik);
%Comprobado
a=[2.03414772537285,0.701742151020526,0.250389962645483]
lik=sum(log(mvnpdf(a,mulik',siglik2)))