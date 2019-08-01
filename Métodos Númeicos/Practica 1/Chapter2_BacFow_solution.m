%upper solution back
x=zeros(n,1);
for j=n:-1:1
    if (U(j,j)==0) error('Matrix is singular!'); end;
    x(j)=b(j)/U(j,j);
    b(1:j-1)=b(1:j-1)-U(1:j-1,j)*x(j);
end
%lower solution foward
for i=1:length(b)
x(i)=(b(i)-A(i,1:i-1)*x(1:i-1))/A(i,i);
end