%Ejercicio 1
%la función a evaluar
fplot(@(z) normal(z), [-3 3])
%function Q=cotest(f,a,b,N,nodes)
%probabilidad trapezoid rule
pr_t=.5+cotest(@normal,0,1,4,4)
fprintf('La probabilidad de que una variable aleatoria sea menor a uno es de: %6.4f\n', pr_t*100)
%probabilidad simpson´s rule
pr_s=.5+cotess(@normal,0,1,4,4);
fprintf('La probabilidad de que una variable aleatoria sea menor a uno es de: %6.4f\n', pr_s*100)
%Dado a que se toma soloun segmento de la distribución normal, se concluye
%que la regla del trapezoide puede ser mejor que con la regla de simpson
b=1;
a=0;
n=4;
h=(b-a)/(n);
w=[h/2;h;h;h/2];
for i=1:n
    x(i)=a+(i-1)*h;
    fx(i)=(1/sqrt(2*pi))*(exp(-x(i)/2));
    fdraw(i)=fx(i)*w(i);
end 
fprintf('res %6.4f\n', sum(fdraw)+0.5) 
   
alpha=1;
alpha_f(1,4)


%Ejercicio 2
%metodo de newton
tic
tol=1e-15;
for alpha=0:.00001:2 
integral=cotest(@alpha_f,0,1,4,4);
 if abs(integral-1)<=tol, break, end
    toc
end 
   

