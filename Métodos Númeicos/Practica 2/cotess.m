%Newton-Cotes function %  Simpson's Rule
function Q=cotess(f,a,b,N,nodes)
%permite hacer la suma de los intervalos que se evaluaran
N=(nodes-1)*ceil(N/(nodes-1));  N1=N+1;
%linea que distribuye los intervalos del mismo tama�o
x=linspace(a,b,N1)';
%largo de cada intervalo
h=x(2)-x(1); 
%funci�n a evaluar
g=f(x);
%puntos finales de la distribuci�na evaluar en la integral
pts_extremos=g(1)+g(N1);
%para determinar la aproximaci�n de la integral mediante la regla de
%Simpson's
%pts_medios son los extremos, y el resto los valores medios
Q=(h/3)*(pts_extremos+4*sum(g(2:2:N))+2*sum(g(3:2:N)));
end
