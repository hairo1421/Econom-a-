%Ej.1
clear all
syms x
%Funciòn exponencial
fplot(1-exp(2*x),[-1,1],':','LineWidth',0.1)
hold on
%Serie de taylor primer orden
fplot(taylor(1-exp(2*x),x,'ExpansionPoint', 0,'order',1),[-1,1],':','LineWidth',0.01)
hold on
%Serie de taylor segundo orden
fplot(taylor(1-exp(2*x),x,'ExpansionPoint', 0,'order',2),[-1,1],':','LineWidth',0.01)
hold on
%Serie de taylor tercer orden
fplot(taylor(1-exp(2*x),x,'ExpansionPoint', 0,'order',3),[-1,1],':','LineWidth',0.01)
hold off
grid on

%Ej. 2 (1.2)
clear all
%a)Normal
A=[0,-1,2;-2,-1,4;2,7,-3];
B=[-7,1,1;7,-3,-2;3,5,0];
Y=[3,-1,2]';
C=A*B;
X=inv(C'*C)*(C'*Y);
%b)Element by Element
A2=[0,-1,2;-2,-1,4;2,7,-3];
B2=[-7,1,1;7,-3,-2;3,5,0];
Y2=[3,-1,2]';
C2=A2.*B2;
X2=inv(C2'*C2)*(C2'*Y2);

%Ej. 3 (1.3)
clear all
t=1960:1:2001;
c=5;
b1=0.05;
e=randn(1,42)*.02;
%Generando la serie
y=c+b1*t+e;
plot(y);
%Regression
kte=ones(1,42);
X=[kte',t'];
b=inv(X'*X)*(X'*y');
%Estimación
y_est=X*b;
plot(t,y)
hold on
plot(t,y_est)
hold off
grid on

%Ej. 4 (2.2)
clear all
A=[54,14,-11,2;14,50,-4,29;-11,-4,55,22;2,29,22,95];
b=[1;1;1;1];
%LU
[L,U] = lu(A)
Y=inv(L)*b
x=inv(U)*Y
%Cholesky
X=A\b;
%Inversa de la matriz A
X2=inv(A)*b;
%% Jacobi Method
maxit=100;
tol=1e-5;
d = diag(A);
x_jm=[0 0 0 0]';
for it=1:maxit
dx = (b-A*x_jm)./d;
x_jm = x_jm+dx;
if norm(dx)<tol, break, end
end
%% Seidel Method
Q = tril(A);
lambda=.5;
x_sm=[0 0 0 0]';
for it=1:maxit
dx2 = Q\(b-A*x_sm);
x_sm = x_sm+lambda*dx2;
if norm(dx2)<tol, break, end
end

