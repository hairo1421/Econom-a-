clear
%DGP AR(2)
T=100;
c=1.5;
b1=.5;
b2=.3;
sigerr=1;
%State space
beta=[b1,b2;1,0];
alpha=[c;0];
%Momentos
mu=inv(eye(2)-beta)*alpha;
sigma=[sigerr,0;0,0];
vecvarcov=inv(eye(2^2)-kron(beta,beta))*sigma(:);
vec=reshape(vecvarcov,2,2);
%creando observaciones
%Nota:Debes ponerlo transpuesto
y=mvnrnd(mu,vec)';
%Generando la serie
for t=3:T
    y(t)=c+b1*y(t-1)+b2*y(t-2)+rand;
end 
 plot(y);
%comprpobra
X=[ones(T-2,1),y(2:end-1),y(1:end-2)];
Y=y(3:end);
ols=inv(X'*X)*X'*Y
 %MCMC
 
 %Momentos Prior
 mu0=[.5;.5;.5];
 sig0=[12,0,0;0,.1,0;0,0,.1^2];
 A0=1/2
 B0=2/2
 %Momentos likelihood
mulik=ols;
siglik=3;

for j=1:8000;
 sigpos=inv(inv(sig0)+inv(siglik));
 mupos=sigpos*(inv(sig0)*mu0+inv(siglik)*mulik);
 theta(j,:)=mvnrnd(mupos,sigpos);
 ssr=sum((Y-X*ols).^2);
 A=A0+ssr/2;
 B=B0+T/2;
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
title('par autoregresivo 1 posterior');
subplot(2,2,4);
%nota:primero filas luego columnas
ksdensity(theta(n0+1:end,3));
title('par autoregresivo 2 posterior');

%HPDI
intsup=mupos(2)+1.96*sqrt(sigpos(5))
intinf=mupos(2)-1.96*sqrt(sigpos(5))
beta1=[intinf,mupos(2),intsup]
intsup=mupos(3)+1.96*sqrt(sigpos(9))
intinf=mupos(3)-1.96*sqrt(sigpos(9))
beta2=[intinf,mupos(3),intsup]
intsup_c=mupos(1)+1.96*sqrt(sigpos(1))
intinf_c=mupos(1)-1.96*sqrt(sigpos(1))
constante=[intinf_c,mupos(1),intsup_c]

%Problemas de Computo
prior=mvnpdf(theta,mu0',sig0)
plot(prior)
title('prior')
prior2=gampdf(s2draw,A0,B0)
siglik2=[siglik,0,0;0,3,0;0,0,3]
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