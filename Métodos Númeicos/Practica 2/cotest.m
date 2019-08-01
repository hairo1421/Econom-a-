%Newton-Cotes function % Trapezoidal Rule
function Q=cotest(f,a,b,N,nodes)
%permite hacer la suma de los intervalos que se evaluaran
N=(nodes-1)*ceil(N/(nodes-1));  N1=N+1;
%linea que distribuye los intervalos del mismo tamaño
x=linspace(a,b,N1)';
%largo de cada intervalo
h=x(2)-x(1);  
%función a evaluar
g=f(x);
%puntos finales de la distribucióna evaluar en la integral
pts_extremos=g(1)+g(N1);
%para determinar la aproximación de la integral mediante la regla del
%trapezoide
%pts_medios son los extremos, y el resto los valores medios
Q=(h/2)*(pts_extremos+2*sum(g(2:N)));
end 
